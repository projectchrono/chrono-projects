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
// This test monitors performance of a single, laminated ANCF humvee tire
// =============================================================================

#include "chrono/ChConfig.h"
#include "chrono/assets/ChVisualShapeCylinder.h"
#include "chrono/assets/ChVisualShapeTriangleMesh.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsValidation.h"
#include "chrono/utils/ChOpenMP.h"

#include "chrono/fea/ChElementShellANCF_3423.h"
#include "chrono/fea/ChLinkNodeSlopeFrame.h"
#include "chrono/fea/ChLinkNodeFrame.h"
#include "chrono/fea/ChMesh.h"

#ifdef CHRONO_PARDISO_MKL
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"
#define USE_MKL
#else
#undef USE_MKL
#endif

#ifdef CHRONO_OPENMP_ENABLED
#include <omp.h>
#endif

using namespace chrono;
using namespace fea;

int num_threads = 4;
double step_size = 0.0002;
int num_steps = 10;

bool addConstRim = true;
bool addBodies = true;
bool addGroundForces = false;
bool addSingleLoad = false;
bool addPressureAlessandro = true;

std::shared_ptr<ChBody> BGround;
std::shared_ptr<ChBody> SimpChassis;           // Chassis body
std::shared_ptr<ChLinkNodeFrame> constraint;  // Create shared pointers for rim-mesh constraints
std::shared_ptr<ChLinkNodeSlopeFrame> constraintD;
std::shared_ptr<ChNodeFEAxyzD> ConstrainedNode;
std::shared_ptr<ChLinkLockPlanar> constraintRim;

auto MloadcontainerGround = chrono_types::make_shared<ChLoadContainer>();
// Some model parameters
const double spring_coef = 3e4;  // Springs and dampers for strut
const double damping_coef = 1e3;
const double ForVelocity = 2;  // Initial velocity of the tire. Applied to hub, nodes, and chassis
double TirePressure =
    -200e3;  // Applied suddently at the beginning of simulation// Was 120e3 // I CHANGED THE SIGN 1/25/2016 9:00PM
const double HumveeVertPos = 0.4673;  // Vertical (z axis) center of rim
// Location of the wheels w.r.t. rigid body's center of mass
const double Lwx = 1.60;     // Value of longitudinal distance
const double Lwy = 1.0;      // Value of lateral distance
const int NoTires = 4;       // Number of tires. Programmer
double GroundLoc = -0.0001;  // Ensure there is no contact with ground during test
double BumpRadius = 0.05;    // changed from 0.035 to 0.05
double BumpLongLoc = 2.3;
 
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
        fscanf(inputfile, "%d %d %d %d %d %d %d\n", &count, &dummy, &SectionID(i, 0), &NodesPerElement(i, 0),
               &NodesPerElement(i, 1), &NodesPerElement(i, 3), &NodesPerElement(i, 2));
        // printf("SectionID[i] %d\n  ", SectionID(i, 0));

        fscanf(inputfile, " %lf %lf\n", &ElementLength(i, 0), &ElementLength(i, 1));
        if (MaxSectionNumber < SectionID(i, 0)) {
            MaxSectionNumber = SectionID(i, 0);
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

    /// "Virtual" copy constructor (covariant return type).
    virtual MyLoadCustomMultiple* Clone() const override { return new MyLoadCustomMultiple(*this); }

    double GroundLocationBump(double GroundLoc, double BumpLoc, ChVector3d NodeLocation, double Amplitude) {
        if (NodeLocation.y() > 0.0 || NodeLocation.x() <= (BumpLoc - BumpRadius) ||
            NodeLocation.x() >= (BumpLoc + BumpRadius))  // There is no bump on that side
        {
            return GroundLoc;
        } else {
            return (GroundLoc +
                    Amplitude * sin(1 / (2 * BumpRadius) * CH_PI * (NodeLocation.x() - (BumpLongLoc - BumpRadius))));
        }
    };
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
        double FrictionCoeff = 0.9;

        // Calculation of number of nodes in contact with the Ground
        int NoCNodes = 0;
        for (int iii = 0; iii < loadables.size(); iii++) {
            Node1_Pos = state_x->segment(iii * 6, 3);
            if (Node1_Pos.z() < GroundLocationBump(GroundLoc, BumpLongLoc, Node1_Pos, BumpRadius)) {
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
                double GroundLocZ = GroundLocationBump(GroundLoc, BumpLongLoc, Node1_Pos, BumpRadius);
                if (Node1_Pos.z() < GroundLocZ) {
                    double Penet = abs(Node1_Pos.z() - GroundLocZ);
                    // std::cout << "Node number:  " << iii << ".  "
                    //          << "Penetration:  " << Penet << "\n";
                    NormalForceNode = KGround * Penet;  // +CGround * abs(Node1_Vel.y()*Penet);
                    this->load_Q(iii * 6 + 2) =
                        NormalForceNode - CGround * (Node1_Vel.z()) * abs(Penet);  // Fy (Vertical)
                    // Friction forces
                    const double VelLimit = 0.25;
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
    virtual bool IsStiff() { return false; }
};

void MakeANCFHumveeWheel(ChSystem& my_system,
                         std::shared_ptr<ChMesh>& TireMesh,
                         const ChVector3d rim_center,
                         std::shared_ptr<ChBody>& Hub_1,
                         double TirePressure,
                         double ForVelocity) {
    // Create rim for this mesh
    my_system.AddBody(Hub_1);
    Hub_1->SetFixed(false);
    Hub_1->EnableCollision(false);
    Hub_1->SetMass(10);
    Hub_1->SetInertiaXX(ChVector3d(0.3, 0.3, 0.3));
    Hub_1->SetPos(rim_center);  // Y = -1m
    Hub_1->SetLinVel(ChVector3d(ForVelocity, 0, 0));
    Hub_1->SetAngVelParent(ChVector3d(0, ForVelocity / (HumveeVertPos),
                                      0));  // 0.3 to be substituted by an actual measure of the average radius.

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
    ChVectorN<int, 2880> SectionID;  // Catagorizes which tire section the elements are a part of
    ChMatrixNM<double, 15, 2> LayPROP;   // Thickness and ply angles of the layered elements
    ChMatrixNM<int, 15, 7> MatID;        // Catagorizes the material of each layer
    ChMatrixNM<double, 7, 12> MPROP;     // Material properties
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
                my_system.Add(constraint);

                // Add rotation constraints
                constraintD = chrono_types::make_shared<ChLinkNodeSlopeFrame>();
                constraintD->Initialize(ConstrainedNode, Hub_1);
                constraintD->SetDirectionInAbsoluteCoords(ConstrainedNode->GetSlope1());
                my_system.Add(constraintD);
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
    my_system.Add(TireMesh);
    auto Mloadcontainer = chrono_types::make_shared<ChLoadContainer>();
    // Add constant pressure using ChLoaderPressure (preferred for simple, constant pressure)
    if (addPressureAlessandro) {
        for (int NoElmPre = 0; NoElmPre < TotalNumElements; NoElmPre++) {
            auto face = std::static_pointer_cast<ChElementShellANCF_3423>(TireMesh->GetElement(NoElmPre));
            auto faceloader = chrono_types::make_shared<ChLoaderPressure>(face);
            faceloader->SetPressure(-TirePressure);
            faceloader->SetStiff(false);
            faceloader->SetIntegrationPoints(2);
            auto faceload = chrono_types::make_shared<ChLoad>(faceloader);
            Mloadcontainer->Add(faceload);
        }
    }
    my_system.Add(Mloadcontainer);
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
        my_system.Add(mloadcontainerGround);
    }  // End loop over tires
}

int main(int argc, char* argv[]) {
    // ---------------------------------
    // Set path to Chrono data directory
    // ---------------------------------
    SetChronoDataPath(CHRONO_DATA_DIR);

    // Definition of the model
    ChSystemNSC my_system;
    utils::ChWriterCSV out("\t");
    out.Stream().setf(std::ios::scientific | std::ios::showpos);
    out.Stream().precision(7);

    // Set number of threads
#ifdef CHRONO_OPENMP_ENABLED
    my_system.SetNumThreads(std::min(num_threads, ChOMP::GetNumProcs()));
#else
    std::cout << "No OpenMP\n";
#endif 

    // Body 1: Ground
    BGround = chrono_types::make_shared<ChBody>();
    my_system.AddBody(BGround);
    BGround->SetFixed(true);
    BGround->EnableCollision(false);
    BGround->SetMass(1);
    BGround->SetInertiaXX(ChVector3d(1, 1, 0.2));
    BGround->SetPos(ChVector3d(-2, 0, 0));  // Y = -1m
    ChQuaternion<> rot = QuatFromAngleX(0.0);
    BGround->SetRot(rot);

    // Create hubs and tire meshes for 4 wheels
    auto Hub_1 = chrono_types::make_shared<ChBody>();
    ChVector3d rim_center_1(0.0, 0.0, HumveeVertPos);  //

    // Create tire meshes
    auto TireMesh1 = chrono_types::make_shared<ChMesh>();
    MakeANCFHumveeWheel(my_system, TireMesh1, rim_center_1, Hub_1, TirePressure, ForVelocity);

    auto mmaterial = chrono_types::make_shared<ChContactMaterialNSC>();
    mmaterial->SetFriction(0.4f);
    mmaterial->SetCompliance(0.0000005f);
    mmaterial->SetComplianceT(0.0000005f);
    mmaterial->SetDampingF(0.2f);
    auto SimpChassis = chrono_types::make_shared<ChBodyAuxRef>();  // visualization?
    my_system.AddBody(SimpChassis);
    SimpChassis->SetMass(2000.0);
    // optional, attach a texture for better visualization
    SimpChassis->SetPos(ChVector3d(0, 0, HumveeVertPos));
    SimpChassis->SetLinVel(ChVector3d(ForVelocity, 0, 0));
    SimpChassis->SetFixed(false);

    auto mtrussmesh = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    mtrussmesh->GetMesh()->LoadWavefrontMesh(GetChronoDataFile("vehicle/hmmwv/hmmwv_chassis_simple.obj"));
    // mtrussmesh->GetMesh()->Transform(VNULL, QuatFromAngleAxis(CH_PI_2, VECT_Z) % QuatFromAngleAxis(CH_PI_2, VECT_Y));
    mtrussmesh->GetMesh()->Transform(VNULL, QuatFromAngleX(0));
    SimpChassis->AddVisualShape(mtrussmesh);

    auto Bump = chrono_types::make_shared<ChBody>();
    Bump->SetMass(10);
    Bump->SetFixed(true);
    Bump->SetPos(ChVector3d(BumpLongLoc, -1.0, 0.0));
    my_system.Add(Bump);

    auto cyl_wheel = chrono_types::make_shared<ChVisualShapeCylinder>(BumpRadius, 1.0);
    Bump->AddVisualShape(cyl_wheel, ChFrame<>(VNULL, QuatFromAngleX(CH_PI_2)));

    // Create joints between chassis and hubs
    auto RevTr_1 = chrono_types::make_shared<ChLinkRevoluteTranslational>();
    my_system.AddLink(RevTr_1);
    RevTr_1->Initialize(Hub_1, SimpChassis, true, ChVector3d(0, 0, 0), ChVector3d(0, 1, 0), ChVector3d(0, 0, 0.1),
                        ChVector3d(0, 0, 1), ChVector3d(1, 0, 0), true);

    // Spring and damper for secondary suspension: True position vectors are relative
    auto spring1 = chrono_types::make_shared<ChLinkTSDA>();
    spring1->Initialize(Hub_1, SimpChassis, true, ChVector3d(0, 0, 0), ChVector3d(Lwx, -Lwy, 0));
    spring1->SetSpringCoefficient(spring_coef);
    spring1->SetDampingCoefficient(damping_coef);
    my_system.AddLink(spring1);

    // Create a large cube as a floor.
    auto mrigidBody = chrono_types::make_shared<ChBodyEasyBox>(20, 20, 0.00001, 1000, true, false); // no collision
    my_system.Add(mrigidBody);
    mrigidBody->SetPos(ChVector3d(0, 0, GroundLoc));
    mrigidBody->SetFixed(true);
    my_system.SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));

// Set up solver
#ifdef USE_MKL
    std::cout << "Using PardisoMKL solver\n";
    auto mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>();
    mkl_solver->LockSparsityPattern(true);
    my_system.SetSolver(mkl_solver);
#else
    std::cout << "Using MINRES solver\n";
    auto solver = chrono_types::make_shared<ChSolverMINRES>();
    solver->EnableDiagonalPreconditioner(true);
    solver->SetMaxIterations(100);
    solver->SetVerbose(false);
    solver->SetTolerance(1e-12);
    my_system.SetSolver(solver);
#endif

    my_system.SetTimestepperType(ChTimestepper::Type::HHT);
    // my_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);  // fast, less precise
    auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(my_system.GetTimestepper());
    mystepper->SetAlpha(-0.2);  // Important for convergence
    mystepper->SetMaxIters(16);
    mystepper->SetAbsTolerances(6e-2, 0.8);
    mystepper->SetModifiedNewton(false);
    mystepper->SetVerbose(true);
    mystepper->SetRequiredSuccessfulSteps(2);
    mystepper->SetMaxItersSuccess(7);

    // Initialize total number of iterations and timer.
    int num_iterations = 0;
    ChTimer timer;
    timer.start();

    // Simulate to final time, while accumulating number of iterations.
    for (int istep = 0; istep < num_steps; istep++) {
        my_system.DoStepDynamics(step_size);
        num_iterations += mystepper->GetNumIterations();
    }
    timer.stop();

    // Report run time and total number of iterations.
    std::cout << "Number of iterations: " << num_iterations << "\n";
    std::cout << "Simulation time:  " << timer() << "\n";
    std::cout << "Internal forces (" << TireMesh1->GetNumCallsInternalForces()
             << "):  " << TireMesh1->GetTimeInternalForces() << "\n";
    std::cout << "Jacobian (" << TireMesh1->GetNumCallsJacobianLoad() << "):  " << TireMesh1->GetTimeJacobianLoad()
             << "\n";
    std::cout << "Extra time:  " << timer() - TireMesh1->GetTimeInternalForces() - TireMesh1->GetTimeJacobianLoad()
             << "\n";
    getchar();

    return 0;
}
