// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Antonio Recuero
// =============================================================================
//
// This test creates four wheels objects using the function makeANCFwheel. The 4
// wheels are constrained to the rim, which in turn is linked to the chassis
// through a ChLinkRevoluteTranslational joint. Values for the spring and damping
// coefficients of the secondary suspension may be selected in the parameter
// definition section.
// =============================================================================

#include "chrono/assets/ChTexture.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsValidation.h"

#include "chrono/fea/ChElementShellANCF_3423.h"
#include "chrono/fea/ChLinkNodeSlopeFrame.h"
#include "chrono/fea/ChLinkNodeFrame.h"
#include "chrono/fea/ChMesh.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"

using namespace chrono;
using namespace fea;
using namespace chrono::irrlicht;

bool addConstRim = true;  //
bool addBodies = true;
bool addGroundForces = true;
bool showVisual = true;
bool addSingleLoad = false;
bool addPressureAlessandro = true;
std::shared_ptr<ChBody> BGround;
std::shared_ptr<ChBodyEasyBox> SimpChassis;    // Chassis body
std::shared_ptr<ChLinkNodeFrame> constraint;  // Create shared pointers for rim-mesh constraints
std::shared_ptr<ChLinkNodeSlopeFrame> constraintD;
std::shared_ptr<ChNodeFEAxyzD> ConstrainedNode;
std::shared_ptr<ChLinkLockPlanar> constraintRim;

auto MloadcontainerGround = chrono_types::make_shared<ChLoadContainer>();
// Some model parameters
const double spring_coef = 3e4;  // Springs and dampers for strut
const double damping_coef = 1e3;
const int num_steps = 1500;    // Number of time steps for unit test (range 1 to 4000) 550
double time_step = 0.00015;    //
const double ForVelocity = 2;  // Initial velocity of the tire. Applied to hub, nodes, and chassis
double TirePressure =
    -120e3;  // Applied suddently at the beginning of simulation// Was 120e3 // I CHANGED THE SIGN 1/25/2016 9:00PM
const double HumveeVertPos = 0.4673;  // Vertical (z axis) center of rim
// Location of the wheels w.r.t. rigid body's center of mass
const double Lwx = 1.1;           // Value of longitudinal distance
const double Lwy = 0.7;           // Value of lateral distance
const int NoTires = 4;            // Number of tires. Programmer
const double GroundLoc = 0.0000;  // Redefined?

// Read input file for the comprehensive Humvee tire
void ReadInputFile(ChMatrixDynamic<double>& COORDFlex,
                   ChMatrixDynamic<double>& VELCYFlex,
                   ChMatrixDynamic<int>& NodesPerElement,
                   int& TotalNumElements,
                   int& NumElements_x,
                   int& NumElements_y,
                   int& TotalNumNodes,
                   ChVectorN<int, 2880>& SectionID,
                   ChMatrixNM<double, 15, 2>& LayerPROP,
                   ChMatrixNM<int, 15, 7>& MatID,
                   ChMatrixNM<double, 7, 12>& MPROP,
                   ChMatrixDynamic<double>& ElementLength,
                   ChVectorN<int, 3>& NumLayPerSect) {
    FILE* inputfile;
    char str1[100];
    int numFlexBody = 0;
    int dummy;
    int count;

    // double ACCELFlex[4000][6];
    // double ACCELRigid[2][7];
    // double LuGreZStart[25][40];
    // double LuGreZStart_dt[25][40];
    // double LuGreZStart_dtdt[25][40];

    int NDR[4000][6];
    int NumLayer[10];
    int MaxSectionNumber = 0;
    int MaxMatID = 0;
    int MTYPE = 0;

    int MAXCOUNT = 100;
    inputfile = fopen(GetChronoDataFile("fea/ANCFtire/IndataBiLinearShell_Tire(HMMWV50x24).INP").c_str(), "r");
    printf("Open IndataBiLinearShell_Tire(HMMWV50x24).INP \n");
    if (inputfile == NULL) {
        printf("Input data file not found!!\n");
        exit(1);
    }

    TotalNumElements = 0;
    NumElements_x = 0;
    NumElements_y = 0;
    TotalNumNodes = 0;

    //!--------------------------------------!
    //!-- Elememt data            -----------!
    //!--------------------------------------!

    fgets(str1, MAXCOUNT, inputfile);
    printf("%s\n", str1);
    fscanf(inputfile, "%d\n", &numFlexBody);

    fgets(str1, MAXCOUNT, inputfile);
    printf("%s\n", str1);
    fscanf(inputfile, "%d %d %d %d\n", &TotalNumElements, &NumElements_x, &NumElements_y, &TotalNumNodes);
    fgets(str1, MAXCOUNT, inputfile);

    printf("%s\n", str1);
    for (int i = 0; i < TotalNumElements; i++) {
        fscanf(inputfile, "%d %d %d %d %d %d %d\n", &count, &dummy, &SectionID(i), &NodesPerElement(i, 0),
               &NodesPerElement(i, 1), &NodesPerElement(i, 3), &NodesPerElement(i, 2));
        printf("SectionID[i] %d\n  ", SectionID(i));

        fscanf(inputfile, " %lf %lf\n", &ElementLength(i, 0), &ElementLength(i, 1));
        if (MaxSectionNumber < SectionID(i)) {
            MaxSectionNumber = SectionID(i);
        }
        // if(TotalNumNodes<max(NumNodes[i][0],max(NumNodes[i][1],max(NumNodes[i][2],NumNodes[i][3]))))
        //{TotalNumNodes=max(NumNodes[i][0],max(NumNodes[i][1],max(NumNodes[i][2],NumNodes[i][3])));}

        // printf("MaxSectionNumber %lf, %lf \n ", ElemLengthXY[i][0],ElemLengthXY[i][1]);
    }

    //!--------------------------------------!
    //!-- NDR,COORDFlex,VELCYFlex -----------!
    //!--------------------------------------!
    // fscanf(inputfile,"%s\n",str1);
    fgets(str1, MAXCOUNT, inputfile);
    printf("%s\n", str1);
    for (int i = 0; i < TotalNumNodes; i++) {
        fscanf(inputfile, "%d %d %d %d %d %d %d\n", &count, &NDR[i][0], &NDR[i][1], &NDR[i][2], &NDR[i][3], &NDR[i][4],
               &NDR[i][5]);
        fscanf(inputfile, "%lf %lf %lf %lf %lf %lf\n", &COORDFlex(i, 0), &COORDFlex(i, 1), &COORDFlex(i, 2),
               &COORDFlex(i, 3), &COORDFlex(i, 4), &COORDFlex(i, 5));
        fscanf(inputfile, "%lf %lf %lf %lf %lf %lf\n", &VELCYFlex(i, 0), &VELCYFlex(i, 1), &VELCYFlex(i, 2),
               &VELCYFlex(i, 3), &VELCYFlex(i, 4), &VELCYFlex(i, 5));
        // printf("NumNodes %d %d %d %d %d %d\n",NDR[i][0],NDR[i][1],NDR[i][2],NDR[i][3],NDR[i][4],NDR[i][5]);
        // printf("NumNodes %lf %lf %lf %lf %lf
        // %lf\n",COORDFlex[i][0],COORDFlex[i][1],COORDFlex[i][2],COORDFlex[i][3],COORDFlex[i][4],COORDFlex[i][5]);
    }

    //!--------------------------------------!
    //!--- Read Layer Data ------------------!
    //!--------------------------------------!
    // fscanf(inputfile,"%s\n",str1);
    fgets(str1, MAXCOUNT, inputfile);
    printf("%s\n", str1);
    int counted = 0;
    for (int i = 0; i < MaxSectionNumber; i++) {
        fscanf(inputfile, "%d %d\n", &count, &NumLayer[i]);
        for (int j = 0; j < NumLayer[i]; j++) {
            fscanf(inputfile, "%lf %lf %d\n", &LayerPROP(counted + j, 0), &LayerPROP(counted + j, 1), &MatID(i, j));
            if (MaxMatID < MatID(i, j)) {
                MaxMatID = MatID(i, j);
            }
            NumLayPerSect(i) = NumLayer[i];
            // printf("%lf %lf %d\n%d\n", LayerPROP(counted + j, 0), LayerPROP(counted + j, 1), MatID(i, j), counted +
            // j);
        }
        counted += NumLayPerSect(i);
    }

    //!--------------------------------------!
    //!--- Read Material Data ---------------!
    //!--------------------------------------!
    // fscanf(inputfile,"%s\n",str1);
    fgets(str1, MAXCOUNT, inputfile);
    printf("%s\n", str1);
    for (int i = 0; i < MaxMatID; i++) {
        fscanf(inputfile, "%d %d\n", &count, &MTYPE);
        if (MTYPE == 1) {
            fscanf(inputfile, "%lf %lf %lf %lf\n", &MPROP(i, 0), &MPROP(i, 1), &MPROP(i, 2), &MPROP(i, 3));
        }
        if (MTYPE == 2) {
            fscanf(inputfile, "%lf %lf %lf %lf\n", &MPROP(i, 0), &MPROP(i, 1), &MPROP(i, 2), &MPROP(i, 3));
            fscanf(inputfile, "%lf %lf %lf %lf %lf %lf\n", &MPROP(i, 4), &MPROP(i, 5), &MPROP(i, 6), &MPROP(i, 7),
                   &MPROP(i, 8), &MPROP(i, 9));
        }
        // printf("%lf %lf %lf %lf\n",MPROP[i][0],MPROP[i][1],MPROP[i][2],MPROP[i][3]);
    }
};

// ChLoadCustomMultiple to include basic node-Ground contact interaction
class MyLoadCustomMultiple : public ChLoadCustomMultiple {
  public:
    MyLoadCustomMultiple(std::vector<std::shared_ptr<ChLoadable>>& mloadables) : ChLoadCustomMultiple(mloadables){};

    virtual MyLoadCustomMultiple* Clone() const override { return new MyLoadCustomMultiple(*this); }

    virtual void ComputeQ(ChState* state_x,      ///< state position to evaluate Q
                          ChStateDelta* state_w  ///< state speed to evaluate Q
    ) {
        std::vector<std::shared_ptr<ChLoadable>> NodeList;

        ChVector3d Node1_Pos;
        ChVector3d Node1_Vel;
        ChVector3d Node1_Grad;
        ChVector3d Node1_GradVel;
        // this->load_Q.FillElem(0);
        double KGround = 9e5;
        double CGround = KGround;
        double NormalForceNode = 0;
        double FrictionCoeff = 0.7;

        // Calculation of number of nodes in contact with the Ground
        int NoCNodes = 0;
        for (int iii = 0; iii < loadables.size(); iii++) {
            Node1_Pos = state_x->segment(iii * 6, 3);
            if (Node1_Pos.z() < GroundLoc) {
                //  std::cout << " \n Node1_Pos.z(): " << Node1_Pos.z() << "\n GroundLoc: " << GroundLoc << "
                //  Number: " << iii;
                NoCNodes++;
            }
        }
        if (NoCNodes > 0) {
            KGround = 9e5 / double(NoCNodes);
            CGround = 0.001 * KGround;
        }
        // std::cout << "  \n"
        //                  << "Nodes into contact:   " << NoCNodes << " \n";
        if (state_x && state_w) {
            for (int iii = 0; iii < loadables.size(); iii++) {
                Node1_Pos = state_x->segment(iii * 6, 3);
                Node1_Grad = state_x->segment(iii * 6 + 3, 3);
                Node1_Vel = state_w->segment(iii * 6, 3);
                Node1_GradVel = state_w->segment(iii * 6 + 3, 3);
                if (Node1_Pos.z() < GroundLoc) {
                    double Penet = abs(Node1_Pos.z() - GroundLoc);
                    // std::cout << "Node number:  " << iii << ".  "
                    //          << "Penetration:  " << Penet << "\n";
                    NormalForceNode = KGround * Penet;  // +CGround * abs(Node1_Vel.y()*Penet);
                    this->load_Q(iii * 6 + 2) =
                        NormalForceNode - CGround * (Node1_Vel.z()) * abs(Penet);  // Fy (Vertical)
                    // Friction forces
                    const double VelLimit = 0.1;
                    if (abs(Node1_Vel.x()) > VelLimit) {
                        this->load_Q(iii * 6 + 0) =
                            -NormalForceNode * FrictionCoeff *
                            (Node1_Vel.x() / sqrt((pow(Node1_Vel.x(), 2) + pow(Node1_Vel.y(), 2))));  // Fx (Plane x)
                    } else {
                        this->load_Q(iii * 6 + 0) =
                            -NormalForceNode * FrictionCoeff * sin(abs(Node1_Vel.x()) * CH_PI_2 / VelLimit) *
                            (Node1_Vel.x() / sqrt((pow(Node1_Vel.x(), 2) + pow(Node1_Vel.y(), 2))));  // Fx (Plane x)
                    }
                    if (abs(Node1_Vel.y()) > VelLimit) {
                        this->load_Q(iii * 6 + 1) =
                            -NormalForceNode * FrictionCoeff *
                            (Node1_Vel.y() / sqrt((pow(Node1_Vel.x(), 2) + pow(Node1_Vel.y(), 2))));  // Fz (Plane y)
                    } else {
                        this->load_Q(iii * 6 + 1) =
                            -NormalForceNode * FrictionCoeff * sin(abs(Node1_Vel.z()) * CH_PI_2 / VelLimit) *
                            (Node1_Vel.y() / sqrt((pow(Node1_Vel.x(), 2) + pow(Node1_Vel.y(), 2))));  // Fz (Plane y)
                    }
                }
            }
        } else {
            // explicit integrators might call ComputeQ(0,0), null pointers mean
            // that we assume current state, without passing state_x for efficiency
            std::cout << "\n This should never happen \n";
        }
    }
    virtual bool IsStiff() { return true; }
};

void MakeANCFHumveeWheel(ChSystem& sys,
                         const ChVector3d rim_center,
                         std::shared_ptr<ChBody>& Hub_1,
                         double TirePressure,
                         double ForVelocity,
                         int Ident) {
    // Create rim for this mesh
    sys.AddBody(Hub_1);
    Hub_1->SetIdentifier(Ident);
    Hub_1->SetFixed(false);
    Hub_1->EnableCollision(false);
    Hub_1->SetMass(10);
    Hub_1->SetInertiaXX(ChVector3d(0.3, 0.3, 0.3));
    Hub_1->SetPos(rim_center);  // Y = -1m
    Hub_1->SetLinVel(ChVector3d(ForVelocity, 0, 0));
    Hub_1->SetAngVelParent(ChVector3d(0, ForVelocity / (HumveeVertPos),
                                      0));  // 0.3 to be substituted by an actual measure of the average radius.

    // Create tire mesh
    auto TireMesh = chrono_types::make_shared<ChMesh>();

    //  Fixing constraints, initial coordinates and velocities
    // READ INPUT DATA AND CREATE ARRAYS

    // Creating arrays for inputting data
    std::cout << "\n-------------------------------------------------\n";
    std::cout << "TEST: ANCF Tire (Fixed),  implicit integration \n\n";

    // Boolean variables to determine which output files are written
    bool output = true;
    bool output1 = true;
    bool output2 = true;
    bool output3 = true;
    bool output4 = true;

    int TotalNumNodes;
    // Matricies to hold the state informatino for the nodes and rigid bodies
    ChMatrixDynamic<double> COORDFlex(3000, 6);
    ChMatrixDynamic<double> VELCYFlex(3000, 6);
    ChMatrixDynamic<double> ACCELFlex(3000, 6);

    ChMatrixNM<double, 2, 7> COORDRigid;
    ChMatrixNM<double, 2, 7> VELCYRigid;
    ChMatrixNM<double, 2, 7> ACCELRigid;

    ChMatrixDynamic<int> NodesPerElement(2880, 4);  // Defines the connectivity between the elements and nodes
    ChMatrixDynamic<double> ElemLength(2880, 2);    // X and Y dimensions of the shell elements

    int TotalNumElements;
    int NumElements_x;
    int NumElements_y;
    ChVectorN<int, 2880> SectionID;     // Catagorizes which tire section the elements are a part of
    ChMatrixNM<double, 15, 2> LayPROP;  // Thickness and ply angles of the layered elements
    ChMatrixNM<int, 15, 7> MatID;       // Catagorizes the material of each layer
    ChMatrixNM<double, 7, 12> MPROP;    // Material properties
    ChVectorN<int, 3> NumLayPerSection;
    double ContactZ = 0.0;  // Vertical location of the flat ground
    ChVector3d NetContact;  // Net contact forces

    // End of declaration of arrays for inputting data

    // Read actual data file
    ReadInputFile(COORDFlex, VELCYFlex, NodesPerElement, TotalNumElements, NumElements_x, NumElements_y, TotalNumNodes,
                  SectionID, LayPROP, MatID, MPROP, ElemLength, NumLayPerSection);
    ///////////////////////////////////////////////////////////////////////////
    // Assign the humvee mesh properties to our ChMesh
    //// Material List (for HMMWV)
    //// i=0: Carcass
    //// i=1: Steel belt in rubber matrix
    //// i=2: Rubber
    ///////////////////////////////////////////////////////////////////////////

    std::vector<std::shared_ptr<ChMaterialShellANCF>> MaterialList(MPROP.rows());
    for (int i = 0; i < MPROP.rows(); i++) {
        double rho = MPROP(i, 0);
        ChVector3d E(MPROP(i, 1), MPROP(i, 2), MPROP(i, 3));
        ChVector3d nu(MPROP(i, 4), MPROP(i, 5), MPROP(i, 6));
        ChVector3d G(MPROP(i, 7), MPROP(i, 8), MPROP(i, 9));
        MaterialList[i] = chrono_types::make_shared<ChMaterialShellANCF>(rho, E, nu, G);
    }

    // Create a set of nodes for the tire based on the input data
    for (int i = 0; i < TotalNumNodes; i++) {
        auto node = chrono_types::make_shared<ChNodeFEAxyzD>(
            ChVector3d(COORDFlex(i, 0) + rim_center.x(), COORDFlex(i, 1) + rim_center.y(), COORDFlex(i, 2)),
            ChVector3d(COORDFlex(i, 3), COORDFlex(i, 4), COORDFlex(i, 5)));
        node->SetPosDt(ChVector3d(VELCYFlex(i, 0), VELCYFlex(i, 1), VELCYFlex(i, 2)));
        node->SetSlope1Dt(ChVector3d(VELCYFlex(i, 3), VELCYFlex(i, 4), VELCYFlex(i, 5)));
        node->SetPosDt2(ChVector3d(ACCELFlex(i, 0), ACCELFlex(i, 1), ACCELFlex(i, 2)));
        node->SetSlope1Dt2(ChVector3d(ACCELFlex(i, 3), ACCELFlex(i, 4), ACCELFlex(i, 5)));
        node->SetMass(0.0);

        TireMesh->AddNode(node);  // Add nodes to the system
    }
    // Check position of the bottom node
    std::cout << "TotalNumNodes: " << TotalNumNodes << "\n\n";
    auto nodetip = std::dynamic_pointer_cast<ChNodeFEAxyzD>(TireMesh->GetNode((TotalNumElements / 2)));
    std::cout << "X : " << nodetip->GetPos().x() << " Y : " << nodetip->GetPos().y() << " Z : " << nodetip->GetPos().z()
             << "\n\n";
    std::cout << "dX : " << nodetip->GetSlope1().x() << " dY : " << nodetip->GetSlope1().y() << " dZ : " << nodetip->GetSlope1().z()
             << "\n\n";

    int LayerHist = 0;  // Number of layers in the previous tire sections

    // Create all elements of the tire
    for (int i = 0; i < TotalNumElements; i++) {
        auto element = chrono_types::make_shared<ChElementShellANCF_3423>();
        element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzD>(TireMesh->GetNode(NodesPerElement(i, 0) - 1)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(TireMesh->GetNode(NodesPerElement(i, 1) - 1)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(TireMesh->GetNode(NodesPerElement(i, 2) - 1)),
                          std::dynamic_pointer_cast<ChNodeFEAxyzD>(TireMesh->GetNode(NodesPerElement(i, 3) - 1)));
        element->SetDimensions(ElemLength(i, 0), ElemLength(i, 1));

        // Determine the section in which the current element resides
        switch (SectionID(i)) {
            // Bead section
            case 1: {
                LayerHist = 0;
                break;
            }
            // Sidewall section
            case 2: {
                LayerHist = NumLayPerSection(0);
                break;
            }
            // Tread section
            case 3: {
                LayerHist = NumLayPerSection(0) + NumLayPerSection(1);
                break;
            }
        }  // End of switch

        // Give material properties to elements as a construction of layers
        for (int j = 0; j < NumLayPerSection(SectionID(i) - 1); j++) {
            element->AddLayer(LayPROP(LayerHist + j, 0), LayPROP(LayerHist + j, 1) * CH_DEG_TO_RAD,
                              MaterialList[MatID(SectionID(i) - 1, j) - 1]);
            // std::cout << "Thickness: " << LayPROP(LayerHist + j, 0) << "  Ply: " << LayPROP(LayerHist + j, 1) << "
            // Mat: " << MatID(SectionID(i) - 1, j) << "\n";
            // std::cout << "Index: " << LayerHist + j << "   PRev: " << LayerHist << "\n";
        }
        element->SetAlphaDamp(0.01);  // 0.005
        TireMesh->AddElement(element);
    }
    // End of assigning properties to TireMesh (ChMesh)
    // Create constraints for the tire and rim
    // Constrain the flexible tire to the rigid rim body.
    if (addConstRim) {
        for (int i = 0; i < TotalNumNodes; i++) {
            if (i < NumElements_x ||
                i >= TotalNumNodes - NumElements_x) {  // Only constrain the nodes at the ends of the bead section
                ConstrainedNode = std::dynamic_pointer_cast<ChNodeFEAxyzD>(TireMesh->GetNode(i));
                // Add position constraints
                constraint = chrono_types::make_shared<ChLinkNodeFrame>();
                constraint->Initialize(ConstrainedNode, Hub_1);
                sys.Add(constraint);

                // Add rotation constraints
                constraintD = chrono_types::make_shared<ChLinkNodeSlopeFrame>();
                constraintD->Initialize(ConstrainedNode, Hub_1);
                constraintD->SetDirectionInAbsoluteCoords(ConstrainedNode->GetSlope1());
                sys.Add(constraintD);
            }
        }
    }

    // END OF INPUT DATA AND CREATE ARRAYS

    // Add initial velocity to the nodes (for rolling)
    for (unsigned int i = 0; i < TireMesh->GetNumNodes(); ++i) {
        ChVector3d node_pos = std::dynamic_pointer_cast<ChNodeFEAxyzD>(TireMesh->GetNode(i))->GetPos();
        double tang_vel = ForVelocity * (node_pos.z()) / (HumveeVertPos);
        ChVector3d NodeVel(tang_vel, 0, 0.0);
        std::dynamic_pointer_cast<ChNodeFEAxyzD>(TireMesh->GetNode(i))->SetPosDt(NodeVel);
    }

    // Switch off mesh class gravity
    TireMesh->SetAutomaticGravity(false);

    // Add the mesh to the system
    sys.Add(TireMesh);
    auto Mloadcontainer = chrono_types::make_shared<ChLoadContainer>();
    // Add constant pressure using ChLoaderPressure (preferred for simple, constant pressure)
    if (addPressureAlessandro) {
        for (int NoElmPre = 0; NoElmPre < TotalNumElements; NoElmPre++) {
            auto faceload = chrono_types::make_shared<ChLoad<ChLoaderPressure>>(
                std::static_pointer_cast<ChElementShellANCF_3423>(TireMesh->GetElement(NoElmPre)));
            faceload->loader.SetPressure(-TirePressure);
            faceload->loader.SetStiff(false);
            faceload->loader.SetIntegrationPoints(2);
            Mloadcontainer->Add(faceload);
        }
    }
    sys.Add(Mloadcontainer);
    // Constraints for each mesh rim
    auto mloadcontainerGround = chrono_types::make_shared<ChLoadContainer>();

    if (addGroundForces) {
        // Select on which nodes we are going to apply a load

        for (int iNode = 0; iNode < TotalNumNodes; iNode++) {
            std::vector<std::shared_ptr<ChLoadable>> NodeList1;
            auto NodeLoad1 = std::dynamic_pointer_cast<ChNodeFEAxyzD>(TireMesh->GetNode(iNode));
            NodeList1.push_back(NodeLoad1);
            auto Mloadcustommultiple1 = chrono_types::make_shared<MyLoadCustomMultiple>(NodeList1);
            mloadcontainerGround->Add(Mloadcustommultiple1);
        }

    }  // End loop over tires
    sys.Add(mloadcontainerGround);

    if (showVisual) {
        auto mvisualizemeshC = chrono_types::make_shared<ChVisualShapeFEA>(TireMesh);
        mvisualizemeshC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_DOT_POS);
        mvisualizemeshC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
        mvisualizemeshC->SetSymbolsThickness(0.005);
        TireMesh->AddVisualShapeFEA(mvisualizemeshC);

        auto mvisualizemeshwire = chrono_types::make_shared<ChVisualShapeFEA>(TireMesh);
        mvisualizemeshwire->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
        mvisualizemeshwire->SetWireframe(true);
        TireMesh->AddVisualShapeFEA(mvisualizemeshwire);

        auto mvisualizemesh = chrono_types::make_shared<ChVisualShapeFEA>(TireMesh);
        // mvisualizemesh->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_SPEED_NORM);
        mvisualizemesh->SetFEMdataType(ChVisualShapeFEA::DataType::ELEM_STRAIN_VONMISES);
        mvisualizemesh->SetColorscaleMinMax(-0.05, 0.05);
        mvisualizemesh->SetSmoothFaces(true);
        TireMesh->AddVisualShapeFEA(mvisualizemesh);
    }
}

int main(int argc, char* argv[]) {
    // Set path to Chrono data directory
    SetChronoDataPath(CHRONO_DATA_DIR);

    // Definition of the model
    ChSystemNSC sys;

    utils::ChWriterCSV out("\t");
    out.Stream().setf(std::ios::scientific | std::ios::showpos);
    out.Stream().precision(7);
    // Main loop for the definition of 4 meshes

    // Body 1: Ground
    BGround = chrono_types::make_shared<ChBody>();
    sys.AddBody(BGround);
    BGround->SetIdentifier(1);
    BGround->SetFixed(true);
    BGround->EnableCollision(false);
    BGround->SetMass(1);
    BGround->SetInertiaXX(ChVector3d(1, 1, 0.2));
    BGround->SetPos(ChVector3d(-2, 0, 0));  // Y = -1m
    ChQuaternion<> rot = QuatFromAngleX(0.0);
    BGround->SetRot(rot);

    // Create hubs and tire meshes for 4 wheels
    auto Hub_1 = chrono_types::make_shared<ChBody>();
    auto Hub_2 = chrono_types::make_shared<ChBody>();
    auto Hub_3 = chrono_types::make_shared<ChBody>();
    auto Hub_4 = chrono_types::make_shared<ChBody>();

    ChVector3d rim_center_1(Lwx, -Lwy, HumveeVertPos);  //
    ChVector3d rim_center_2(Lwx, Lwy, HumveeVertPos);
    ChVector3d rim_center_3(-Lwx, Lwy, HumveeVertPos);
    ChVector3d rim_center_4(-Lwx, -Lwy, HumveeVertPos);

    MakeANCFHumveeWheel(sys, rim_center_1, Hub_1, TirePressure, ForVelocity, 2);
    MakeANCFHumveeWheel(sys, rim_center_2, Hub_2, TirePressure, ForVelocity, 3);
    MakeANCFHumveeWheel(sys, rim_center_3, Hub_3, TirePressure, ForVelocity, 4);
    MakeANCFHumveeWheel(sys, rim_center_4, Hub_4, TirePressure, ForVelocity, 5);

    auto mmaterial = chrono_types::make_shared<ChContactMaterialNSC>();
    mmaterial->SetFriction(0.4f);
    mmaterial->SetCompliance(0.0000005f);
    mmaterial->SetComplianceT(0.0000005f);
    mmaterial->SetDampingF(0.2f);

    SimpChassis = chrono_types::make_shared<ChBodyEasyBox>(2.4, 1.1, 0.2,  // x,y,z size
                                                           2800,           // density
                                                           true,           // visualization?
                                                           false);         // collision?
    sys.AddBody(SimpChassis);
    // optional, attach a texture for better visualization
    SimpChassis->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/cubetexture_bluewhite.png"));
    auto mtexturebox = chrono_types::make_shared<ChTexture>();
    SimpChassis->SetPos(ChVector3d(0, 0, HumveeVertPos));
    SimpChassis->SetPosDt(ChVector3d(ForVelocity, 0, 0));
    SimpChassis->SetFixed(false);
    // */
    // Create joints between chassis and hubs
    auto RevTr_1 = chrono_types::make_shared<ChLinkRevoluteTranslational>();
    sys.AddLink(RevTr_1);
    RevTr_1->Initialize(Hub_1, SimpChassis, true, ChVector3d(0, 0, 0), ChVector3d(0, 1, 0), ChVector3d(Lwx, -Lwy, 0.1),
                        ChVector3d(0, 0, 1), ChVector3d(1, 0, 0), true);

    auto RevTr_2 = chrono_types::make_shared<ChLinkRevoluteTranslational>();
    sys.AddLink(RevTr_2);
    RevTr_2->Initialize(Hub_2, SimpChassis, true, ChVector3d(0, 0, 0), ChVector3d(0, 1, 0), ChVector3d(Lwx, Lwy, 0.1),
                        ChVector3d(0, 0, 1), ChVector3d(1, 0, 0), true);

    auto RevTr_3 = chrono_types::make_shared<ChLinkRevoluteTranslational>();
    sys.AddLink(RevTr_3);
    RevTr_3->Initialize(Hub_3, SimpChassis, true, ChVector3d(0, 0, 0), ChVector3d(0, 1, 0), ChVector3d(-Lwx, Lwy, 0.1),
                        ChVector3d(0, 0, 1), ChVector3d(1, 0, 0), true);

    auto RevTr_4 = chrono_types::make_shared<ChLinkRevoluteTranslational>();
    sys.AddLink(RevTr_4);
    RevTr_4->Initialize(Hub_4, SimpChassis, true, ChVector3d(0, 0, 0), ChVector3d(0, 1, 0), ChVector3d(-Lwx, -Lwy, 0.1),
                        ChVector3d(0, 0, 1), ChVector3d(1, 0, 0), true);

    // Spring and damper for secondary suspension: True position vectors are relative
    auto spring1 = chrono_types::make_shared<ChLinkTSDA>();
    spring1->Initialize(Hub_1, SimpChassis, true, ChVector3d(0, 0, 0), ChVector3d(Lwx, -Lwy, 0));
    spring1->SetSpringCoefficient(spring_coef);
    spring1->SetDampingCoefficient(damping_coef);
    sys.AddLink(spring1);

    auto spring2 = chrono_types::make_shared<ChLinkTSDA>();
    spring2->Initialize(Hub_2, SimpChassis, true, ChVector3d(0, 0, 0), ChVector3d(-Lwx, -Lwy, 0));
    spring2->SetSpringCoefficient(spring_coef);
    spring2->SetDampingCoefficient(damping_coef);
    sys.AddLink(spring2);

    auto spring3 = chrono_types::make_shared<ChLinkTSDA>();
    spring3->Initialize(Hub_3, SimpChassis, true, ChVector3d(0, 0, 0), ChVector3d(-Lwx, Lwy, 0));
    spring3->SetSpringCoefficient(spring_coef);
    spring3->SetDampingCoefficient(damping_coef);
    sys.AddLink(spring3);

    auto spring4 = chrono_types::make_shared<ChLinkTSDA>();
    spring4->Initialize(Hub_4, SimpChassis, true, ChVector3d(0, 0, 0), ChVector3d(Lwx, Lwy, 0));
    spring4->SetSpringCoefficient(spring_coef);
    spring4->SetDampingCoefficient(damping_coef);
    sys.AddLink(spring4);

    // Create a large cube as a floor.

    auto mrigidBody = chrono_types::make_shared<ChBodyEasyBox>(10, 10, 0.00001, 1000, true, false);
    sys.Add(mrigidBody);
    mrigidBody->SetPos(ChVector3d(0, 0, GroundLoc));
    mrigidBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
    mrigidBody->SetFixed(true);

    sys.SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));
    auto mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>();
    sys.SetSolver(mkl_solver);
    mkl_solver->LockSparsityPattern(true);

    sys.SetTimestepperType(ChTimestepper::Type::HHT);
    // sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);  // fast, less precise
    auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(sys.GetTimestepper());
    mystepper->SetAlpha(-0.3);  // Important for convergence
    mystepper->SetMaxIters(20);
    mystepper->SetAbsTolerances(6e-03, 2.5);
    mystepper->SetModifiedNewton(false);
    mystepper->SetVerbose(true);
    mystepper->SetRequiredSuccessfulSteps(2);
    mystepper->SetMaxItersSuccess(7);

    // Visualization
    /*
    auto mobjmesh = chrono_types::make_shared<ChObjShapeFile>();
    mobjmesh->SetFilename(GetChronoDataFile("fea/tractor_wheel_rim.obj"));
    Hub_1->AddAsset(mobjmesh);
    Hub_2->AddAsset(mobjmesh);
    Hub_3->AddAsset(mobjmesh);
    Hub_4->AddAsset(mobjmesh);
    */

    double start = std::clock();

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(1080, 800);
    vis->SetWindowTitle("HMMWV ANCF tire");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddCamera(ChVector3d(0.5, 0.5, 1.15), ChVector3d(0.65, 0, 0));
    vis->AddTypicalLights();
    vis->AttachSystem(&sys);

    sys.Setup();
    sys.Update();

    std::cout << "\n\nREADME\n\n"
              << " - Press SPACE to start dynamic simulation \n - Press F10 for nonlinear statics - Press F11 for "
                 "linear statics. \n";

    int AccuNoIterations = 0;
    double ChTime = 0.0;

    const double VerForce = 0;
    const double HorForce = 80;
    const double tini = 0.1;
    const double tend = 0.2;
    const double interval = tend - tini;

    while (vis->Run()) {
        vis->BeginScene();
        vis->Render();
        vis->EndScene();

        sys.DoStepDynamics(time_step);

        std::cout << "Time t = " << sys.GetChTime() << "s \n";
        // AccuNoIterations += mystepper->GetNumIterations();
        printf("Vertical position of Tires:      %12.4e       %12.4e       %12.4e       %12.4e  Chassis   \n",
               Hub_1->GetPos().y(), Hub_2->GetPos().y(), Hub_3->GetPos().y(), Hub_4->GetPos().y());

        printf("Longitudinal position of Tires:      %12.4e       %12.4e       %12.4e       %12.4e  Chassis  ",
               Hub_1->GetPos().z(), Hub_2->GetPos().z(), Hub_3->GetPos().z(), Hub_4->GetPos().z());
        out << sys.GetChTime() << Hub_1->GetPos().x() << Hub_1->GetPos().y() << Hub_1->GetPos().z()
            << Hub_2->GetPos().x() << Hub_2->GetPos().y() << Hub_2->GetPos().z() << Hub_3->GetPos().x()
            << Hub_3->GetPos().y() << Hub_3->GetPos().z() << std::endl;
        out.WriteToFile("../VertPosRim.txt");
    }

    double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Computation Time: " << duration;

    /*
    ChVectorDynamic<double> Cp;
    ChVectorDynamic<double> Cd;  // Matrices for storing constraint violations
    sys.Setup();
    sys.Update();
    sys.SetupInitial();
    for (unsigned int it = 0; it < num_steps; it++) {
        sys.DoStepDynamics(time_step);
        std::cout << "Time t = " << sys.GetChTime() << "s \n";
        // std::cout << "nodetip->pos.z() = " << nodetip->pos.z() << "\n";
        // std::cout << "mystepper->GetNumIterations()= " << mystepper->GetNumIterations() << "\n";
        if (addConstRim) {
            Cp = NodePosRim->GetC();
            printf("Point constraint violations:      %12.4e  %12.4e  %12.4e\n", Cp(0), Cp(1), Cp(2));
            Cd = NodeDirRim->GetC();
            printf("Direction constraint violations:  %12.4e  %12.4e\n", Cd(0), Cd(1));

            printf("Vertical position of the rim:  %12.4e m  %12.4e m  %12.4e m  %12.4e m\n", Hub_1->coord.pos.y(),
                   Hub_2->coord.pos.y(), Hub_3->coord.pos.y(), Hub_4->coord.pos.y());
        }
    }
    double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    std::cout << "Computation Time: " << duration;
    */
    return 0;
}
