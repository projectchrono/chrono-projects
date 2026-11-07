// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2026 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Author: Radu Serban
// =============================================================================
//
// Utilities shared by the DEM validation programs (Chrono::Multicore, SMC):
// - contacts of the last step, with forces and stored tangential displacements
// - particle and contact snapshots (CSV)
// - run metadata (JSON)
// - saving and restoring a prepared sample, including its contact history
//
// All quantities are SI. Forces are those acting on the second body of a contact.
// =============================================================================

#ifndef DEM_VALIDATION_UTILS_H
#define DEM_VALIDATION_UTILS_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "chrono/physics/ChBody.h"
#include "chrono_multicore/physics/ChSystemMulticore.h"

namespace dem_validation {

// -----------------------------------------------------------------------------
// Contacts
// -----------------------------------------------------------------------------

/// One contact between two shapes in the last step (only touching contacts are reported).
struct Contact {
    int b1, b2;             // body indices
    int s1, s2;             // global shape indices
    int l1, l2;             // shape indices within the collision model of each body
    chrono::ChVector3d p1;  // contact point on shape 1
    chrono::ChVector3d p2;  // contact point on shape 2
    chrono::ChVector3d n;   // unit normal, from shape 1 towards shape 2
    double overlap;         // penetration (positive)
    double fn;              // normal force on body 2, along n (positive: repulsive)
    chrono::ChVector3d ft;  // tangential force on body 2
    double mu;              // composite sliding friction coefficient
    chrono::ChVector3d ut;  // stored tangential displacement of body 2 relative to body 1 (MultiStep)
    bool history;           // true if a contact history slot holds this contact
};

/// Collect the touching contacts of the last step, with their forces and stored tangential displacements.
inline std::vector<Contact> GetContacts(chrono::ChSystemMulticoreSMC& sys) {
    using namespace chrono;
    std::vector<Contact> contacts;
    auto dm = sys.data_manager;
    auto cd = dm->cd_data;
    if (!cd)
        return contacts;
    const auto& hd = dm->host_data;
    uint nc = cd->num_rigid_contacts;
    if (nc == 0 || hd.ct_force.size() != 2 * (size_t)nc)
        return contacts;

    bool multistep =
        dm->settings.solver.tangential_displ_mode == ChSystemSMC::TangentialDisplacementModel::MultiStep && hd.shear_neigh.size() == (size_t)max_shear * dm->num_rigid_bodies;

    contacts.reserve(nc);
    for (uint i = 0; i < nc; i++) {
        if (cd->dpth_rigid_rigid[i] >= 0)
            continue;
        Contact c;
        c.b1 = cd->bids_rigid_rigid[i].x;
        c.b2 = cd->bids_rigid_rigid[i].y;
        long long sid = cd->contact_shapeIDs[i];
        c.s1 = int(sid >> 32);
        c.s2 = int(sid & 0xffffffff);
        c.l1 = cd->shape_data.local_rigid[c.s1];
        c.l2 = cd->shape_data.local_rigid[c.s2];
        const auto& pa = cd->cpta_rigid_rigid[i];
        const auto& pb = cd->cptb_rigid_rigid[i];
        const auto& nn = cd->norm_rigid_rigid[i];
        c.p1 = ChVector3d(pa.x, pa.y, pa.z);
        c.p2 = ChVector3d(pb.x, pb.y, pb.z);
        c.n = ChVector3d(nn.x, nn.y, nn.z);
        c.overlap = -cd->dpth_rigid_rigid[i];
        const auto& f = hd.ct_force[2 * i + 1];
        ChVector3d F(f.x, f.y, f.z);
        c.fn = F.Dot(c.n);
        c.ft = F - c.fn * c.n;
        c.mu = hd.fric_rigid_rigid[i].x;
        c.ut = VNULL;
        c.history = false;

        // The history of a contact is stored on the body with the larger index, keyed by the other body and the two
        // shapes (larger shape index first). The stored vector is the displacement of shear_body1's partner, with the
        // sign convention of ChIterativeSolverMulticoreSMC.
        if (multistep) {
            int sb1 = std::max(c.b1, c.b2);
            int sb2 = std::min(c.b1, c.b2);
            int ss1 = std::max(c.s1, c.s2);
            int ss2 = std::min(c.s1, c.s2);
            for (int k = 0; k < max_shear; k++) {
                const auto& nb = hd.shear_neigh[max_shear * sb1 + k];
                if (nb.x == sb2 && nb.y == ss1 && nb.z == ss2) {
                    const auto& u = hd.shear_disp[max_shear * sb1 + k];
                    ChVector3d U(u.x, u.y, u.z);
                    c.ut = (sb1 == c.b1) ? U : -U;
                    c.history = true;
                    break;
                }
            }
        }
        contacts.push_back(c);
    }
    return contacts;
}

/// Ratio |f_t| / (mu f_n) of a contact (1 at the Coulomb limit).
inline double CoulombRatio(const Contact& c) {
    return (c.fn > 0 && c.mu > 0) ? c.ft.Length() / (c.mu * c.fn) : 0.0;
}

/// Contact statistics over the particles (bodies with index >= first_particle).
struct ContactStats {
    int num_contacts = 0;     // touching contacts, all
    int num_pp = 0;           // particle-particle contacts
    int num_pw = 0;           // particle-wall contacts
    double coord = 0;         // mean number of contacts per particle (walls included)
    double coord_mech = 0;    // same, over particles with at least two contacts
    int num_rattlers = 0;     // particles with fewer than two contacts
    double frac_sliding = 0;  // fraction of contacts with |f_t| >= 0.99 mu f_n
    double mean_fn = 0;       // mean normal force
};

inline ContactStats ComputeContactStats(const std::vector<Contact>& contacts, int first_particle, int num_particles, std::vector<int>* counts = nullptr) {
    ContactStats st;
    std::vector<int> nc(num_particles, 0);
    int sliding = 0;
    double sum_fn = 0;
    for (const auto& c : contacts) {
        bool p1 = c.b1 >= first_particle;
        bool p2 = c.b2 >= first_particle;
        if (!p1 && !p2)
            continue;
        st.num_contacts++;
        if (p1 && p2)
            st.num_pp++;
        else
            st.num_pw++;
        if (p1)
            nc[c.b1 - first_particle]++;
        if (p2)
            nc[c.b2 - first_particle]++;
        if (CoulombRatio(c) >= 0.99)
            sliding++;
        sum_fn += c.fn;
    }
    long total = 0, total_mech = 0;
    int n_mech = 0;
    for (int k : nc) {
        total += k;
        if (k >= 2) {
            total_mech += k;
            n_mech++;
        } else {
            st.num_rattlers++;
        }
    }
    st.coord = num_particles > 0 ? double(total) / num_particles : 0;
    st.coord_mech = n_mech > 0 ? double(total_mech) / n_mech : 0;
    st.frac_sliding = st.num_contacts > 0 ? double(sliding) / st.num_contacts : 0;
    st.mean_fn = st.num_contacts > 0 ? sum_fn / st.num_contacts : 0;
    if (counts)
        *counts = nc;
    return st;
}

// -----------------------------------------------------------------------------
// Snapshots
// -----------------------------------------------------------------------------

/// Write the state of the particles (one row per particle) to a CSV file, under the given comment line.
inline void WriteParticles(const std::string& filename,
                           const std::vector<std::shared_ptr<chrono::ChBody>>& particles,
                           const std::vector<double>& radii,
                           const std::vector<int>& counts,
                           const std::string& comment) {
    std::ofstream f(filename);
    f << "# " << comment << "\n";
    f << "id,r,x,y,z,vx,vy,vz,wx,wy,wz,q0,q1,q2,q3,nc\n";
    f << std::setprecision(9);
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& b = particles[i];
        const auto& p = b->GetPos();
        const auto& v = b->GetPosDt();
        auto w = b->GetAngVelParent();
        const auto& q = b->GetRot();
        f << b->GetIndex() << "," << radii[i] << "," << p.x() << "," << p.y() << "," << p.z() << "," << v.x() << "," << v.y() << "," << v.z() << "," << w.x() << "," << w.y() << ","
          << w.z() << "," << q.e0() << "," << q.e1() << "," << q.e2() << "," << q.e3() << "," << (counts.empty() ? 0 : counts[i]) << "\n";
    }
}

/// Comment line of the direct shear snapshots: the time and the shear displacement.
inline std::string ShearComment(double time, double shear_disp) {
    std::ostringstream c;
    c << "time = " << std::setprecision(8) << time << " s, shear displacement = " << shear_disp << " m";
    return c.str();
}

/// Write the state of the particles (one row per particle) to a CSV file.
inline void WriteParticles(const std::string& filename,
                           const std::vector<std::shared_ptr<chrono::ChBody>>& particles,
                           const std::vector<double>& radii,
                           const std::vector<int>& counts,
                           double time,
                           double shear_disp) {
    WriteParticles(filename, particles, radii, counts, ShearComment(time, shear_disp));
}

/// Write the touching contacts (one row per contact) to a CSV file, under the given comment line.
inline void WriteContacts(const std::string& filename, const std::vector<Contact>& contacts, const std::string& comment) {
    std::ofstream f(filename);
    f << "# " << comment << "\n";
    f << "# forces act on body b2; n points from shape 1 to shape 2; ut is the stored tangential displacement of b2 "
         "relative to b1 (MultiStep)\n";
    f << "b1,b2,l1,l2,px,py,pz,nx,ny,nz,overlap,fn,ftx,fty,ftz,mu,ratio,utx,uty,utz,hist\n";
    f << std::setprecision(8);
    for (const auto& c : contacts) {
        auto p = 0.5 * (c.p1 + c.p2);
        f << c.b1 << "," << c.b2 << "," << c.l1 << "," << c.l2 << "," << p.x() << "," << p.y() << "," << p.z() << "," << c.n.x() << "," << c.n.y() << "," << c.n.z() << ","
          << c.overlap << "," << c.fn << "," << c.ft.x() << "," << c.ft.y() << "," << c.ft.z() << "," << c.mu << "," << CoulombRatio(c) << "," << c.ut.x() << "," << c.ut.y() << ","
          << c.ut.z() << "," << (c.history ? 1 : 0) << "\n";
    }
}

/// Write the touching contacts (one row per contact) to a CSV file.
inline void WriteContacts(const std::string& filename, const std::vector<Contact>& contacts, double time, double shear_disp) {
    WriteContacts(filename, contacts, ShearComment(time, shear_disp));
}

// -----------------------------------------------------------------------------
// Run metadata
// -----------------------------------------------------------------------------

/// Minimal writer for a flat JSON object, with nested objects passed as JSON text.
class JsonObject {
  public:
    void Add(const std::string& key, const std::string& value) { m_items.push_back({key, Quote(value)}); }
    void Add(const std::string& key, const char* value) { Add(key, std::string(value)); }
    void Add(const std::string& key, bool value) { m_items.push_back({key, value ? "true" : "false"}); }
    void Add(const std::string& key, int value) { m_items.push_back({key, std::to_string(value)}); }
    void Add(const std::string& key, unsigned int value) { m_items.push_back({key, std::to_string(value)}); }
    void Add(const std::string& key, long long value) { m_items.push_back({key, std::to_string(value)}); }
    void Add(const std::string& key, unsigned long long value) { m_items.push_back({key, std::to_string(value)}); }
    void Add(const std::string& key, double value) { m_items.push_back({key, Number(value)}); }
    void Add(const std::string& key, const std::vector<double>& values) {
        std::string s = "[";
        for (size_t i = 0; i < values.size(); i++)
            s += (i ? ", " : "") + Number(values[i]);
        m_items.push_back({key, s + "]"});
    }
    void AddRaw(const std::string& key, const std::string& json) { m_items.push_back({key, json}); }

    std::string Str(int indent = 0) const {
        std::string pad(indent + 2, ' ');
        std::string s = "{\n";
        for (size_t i = 0; i < m_items.size(); i++)
            s += pad + Quote(m_items[i].first) + ": " + m_items[i].second + (i + 1 < m_items.size() ? ",\n" : "\n");
        return s + std::string(indent, ' ') + "}";
    }

    void Write(const std::string& filename) const {
        std::ofstream f(filename);
        f << Str() << "\n";
    }

    static std::string Quote(const std::string& s) {
        std::string q = "\"";
        for (char c : s) {
            if (c == '"' || c == '\\')
                q += '\\';
            q += c;
        }
        return q + "\"";
    }

    static std::string Number(double v) {
        if (!std::isfinite(v))
            return "null";
        std::ostringstream o;
        o << std::setprecision(12) << v;
        return o.str();
    }

  private:
    std::vector<std::pair<std::string, std::string>> m_items;
};

/// Current local date and time as text.
inline std::string Now() {
    std::time_t t = std::time(nullptr);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
    return buf;
}

// -----------------------------------------------------------------------------
// Prepared sample: save and restore, including the Multicore contact history
// -----------------------------------------------------------------------------

/// Save the state of all bodies and the contact history arrays to a binary file.
/// The bodies must be re-created in the same order before the state is restored. 'extra' holds program values.
inline bool SaveState(chrono::ChSystemMulticoreSMC& sys, const std::string& filename, const std::vector<double>& extra) {
    std::ofstream f(filename, std::ios::binary);
    if (!f)
        return false;
    auto wr = [&f](double v) { f.write(reinterpret_cast<const char*>(&v), sizeof(double)); };
    auto wi = [&f](int64_t v) { f.write(reinterpret_cast<const char*>(&v), sizeof(int64_t)); };

    const char magic[8] = {'C', 'H', 'D', 'E', 'M', 'S', 'T', '1'};
    f.write(magic, 8);
    wr(sys.GetChTime());
    const auto& bodies = sys.GetBodies();
    wi((int64_t)bodies.size());
    for (const auto& b : bodies) {
        const auto& p = b->GetPos();
        const auto& q = b->GetRot();
        const auto& v = b->GetPosDt();
        auto w = b->GetAngVelParent();
        for (double x : {p.x(), p.y(), p.z(), q.e0(), q.e1(), q.e2(), q.e3(), v.x(), v.y(), v.z(), w.x(), w.y(), w.z()})
            wr(x);
    }
    const auto& hd = sys.data_manager->host_data;
    wi((int64_t)hd.shear_neigh.size());
    for (size_t i = 0; i < hd.shear_neigh.size(); i++) {
        wi(hd.shear_neigh[i].x);
        wi(hd.shear_neigh[i].y);
        wi(hd.shear_neigh[i].z);
        wr(hd.shear_disp[i].x);
        wr(hd.shear_disp[i].y);
        wr(hd.shear_disp[i].z);
        wr(hd.contact_relvel_init[i]);
        wr(hd.contact_duration[i]);
    }
    wi((int64_t)extra.size());
    for (double x : extra)
        wr(x);
    return (bool)f;
}

/// Restore a state saved with SaveState. Return false if the file does not match the system.
inline bool LoadState(chrono::ChSystemMulticoreSMC& sys, const std::string& filename, std::vector<double>& extra) {
    using namespace chrono;
    std::ifstream f(filename, std::ios::binary);
    if (!f)
        return false;
    auto rd = [&f]() {
        double v = 0;
        f.read(reinterpret_cast<char*>(&v), sizeof(double));
        return v;
    };
    auto ri = [&f]() {
        int64_t v = 0;
        f.read(reinterpret_cast<char*>(&v), sizeof(int64_t));
        return v;
    };

    char magic[8];
    f.read(magic, 8);
    if (std::string(magic, 8) != "CHDEMST1")
        return false;
    double time = rd();
    const auto& bodies = sys.GetBodies();
    if (ri() != (int64_t)bodies.size())
        return false;
    for (const auto& b : bodies) {
        double x[13];
        for (int k = 0; k < 13; k++)
            x[k] = rd();
        b->SetPos(ChVector3d(x[0], x[1], x[2]));
        b->SetRot(ChQuaterniond(x[3], x[4], x[5], x[6]));
        b->SetPosDt(ChVector3d(x[7], x[8], x[9]));
        b->SetAngVelParent(ChVector3d(x[10], x[11], x[12]));
    }
    auto& hd = sys.data_manager->host_data;
    int64_t ns = ri();
    hd.shear_neigh.resize(ns);
    hd.shear_disp.resize(ns);
    hd.contact_relvel_init.resize(ns);
    hd.contact_duration.resize(ns);
    for (int64_t i = 0; i < ns; i++) {
        hd.shear_neigh[i].x = (int)ri();
        hd.shear_neigh[i].y = (int)ri();
        hd.shear_neigh[i].z = (int)ri();
        hd.shear_disp[i].x = rd();
        hd.shear_disp[i].y = rd();
        hd.shear_disp[i].z = rd();
        hd.contact_relvel_init[i] = rd();
        hd.contact_duration[i] = rd();
    }
    int64_t ne = ri();
    extra.resize(ne);
    for (int64_t i = 0; i < ne; i++)
        extra[i] = rd();
    sys.SetChTime(time);
    return (bool)f;
}

}  // end namespace dem_validation

#endif
