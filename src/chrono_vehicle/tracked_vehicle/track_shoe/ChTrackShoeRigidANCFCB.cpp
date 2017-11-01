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
// Authors: Radu Serban, Michael Taylor
// =============================================================================
//
// Base class for a continuous band rigid-link track shoe (template definition).
//
// =============================================================================

#include "chrono/physics/ChGlobal.h"
#include "chrono/assets/ChCylinderShape.h"
#include "chrono/assets/ChBoxShape.h"
#include "chrono/assets/ChColorAsset.h"
#include "chrono/assets/ChTexture.h"

#include "chrono/physics/ChLoadsBody.h"
#include "chrono/physics/ChLoadContainer.h"

#include "chrono_vehicle/ChSubsysDefs.h"
#include "chrono_vehicle/tracked_vehicle/track_shoe/ChTrackShoeRigidANCFCB.h"

#include "chrono_fea/ChElementShellANCF.h"
#include "chrono_fea/ChElementShellANCF_8.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChContactSurfaceMesh.h"
#include "chrono_fea/ChContactSurfaceNodeCloud.h"
#include "chrono_fea/ChNodeFEAbase.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"

#define USE_ANCF_4;
//#define USE_ANCF_8;

using namespace chrono::fea;

namespace chrono {
namespace vehicle {

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

// Utility function to calculate the center of a circle of given radius which
// passes through two given points.
ChVector2<> CalcCircleCenter_RigidANCF(const ChVector2<>& A, const ChVector2<>& B, double r, double direction) {
    // midpoint
    ChVector2<> C = (A + B) / 2;
    // distance between A and B
    double l = (B - A).Length();
    // distance between C and O
    double d = std::sqrt(r * r - l * l / 4);
    // slope of line AB
    double mAB = (B.y() - A.y()) / (B.x() - A.x());
    // slope of line CO (perpendicular to AB)
    double mCO = -1 / mAB;
    // x offset from C
    double x_offset = d / std::sqrt(1 + mCO * mCO);
    // y offset from C
    double y_offset = mCO * x_offset;
    // circle center
    ChVector2<> O(C.x() + direction * x_offset, C.y() + direction * y_offset);

    ////std::cout << std::endl;
    ////std::cout << "radius: " << r << std::endl;
    ////std::cout << A.x() << "  " << A.y() << std::endl;
    ////std::cout << B.x() << "  " << B.y() << std::endl;
    ////std::cout << O.x() << "  " << O.y() << std::endl;
    ////std::cout << "Check: " << (A - O).Length() - r << "  " << (B - O).Length() - r << std::endl;
    ////std::cout << std::endl;

    return O;
}

ChTrackShoeRigidANCFCB::ChTrackShoeRigidANCFCB(const std::string& name) : ChTrackShoe(name) {}

void ChTrackShoeRigidANCFCB::Initialize(std::shared_ptr<ChBodyAuxRef> chassis,
                                    const ChVector<>& location,
                                    const ChQuaternion<>& rotation) {
    // Cache values calculated from template parameters.
    m_seg_length = GetWebLength() / GetNumWebSegments();
    m_seg_mass = GetWebMass() / GetNumWebSegments();
    m_seg_inertia = GetWebInertia();  //// TODO - properly distribute web inertia

    // Cache the postive (+x) tooth arc position and arc starting and ending angles
    ChVector2<> tooth_base_p(GetToothBaseLength() / 2, GetWebThickness() / 2);
    ChVector2<> tooth_tip_p(GetToothTipLength() / 2, GetToothHeight() + GetWebThickness() / 2);
    m_center_p = CalcCircleCenter_RigidANCF(tooth_base_p, tooth_tip_p, GetToothArcRadius(), -1);
    m_center_p_arc_start = std::atan2(tooth_base_p.y() - m_center_p.y(), tooth_base_p.x() - m_center_p.x());
    m_center_p_arc_start = m_center_p_arc_start < 0 ? m_center_p_arc_start + CH_C_2PI : m_center_p_arc_start;
    m_center_p_arc_end = std::atan2(tooth_tip_p.y() - m_center_p.y(), tooth_tip_p.x() - m_center_p.x());
    m_center_p_arc_end = m_center_p_arc_end < 0 ? m_center_p_arc_end + CH_C_2PI : m_center_p_arc_end;
    if (m_center_p_arc_start > m_center_p_arc_end) {
        double temp = m_center_p_arc_start;
        m_center_p_arc_start = m_center_p_arc_end;
        m_center_p_arc_end = temp;
    }

    // Cache the negative (-x) tooth arc position and arc starting and ending angles
    ChVector2<> tooth_base_m(-GetToothBaseLength() / 2, GetWebThickness() / 2);
    ChVector2<> tooth_tip_m(-GetToothTipLength() / 2, GetToothHeight() + GetWebThickness() / 2);
    m_center_m = CalcCircleCenter_RigidANCF(tooth_base_m, tooth_tip_m, GetToothArcRadius(), +1);
    m_center_m_arc_start = std::atan2(tooth_base_m.y() - m_center_m.y(), tooth_base_m.x() - m_center_m.x());
    m_center_m_arc_start = m_center_m_arc_start < 0 ? m_center_m_arc_start + CH_C_2PI : m_center_m_arc_start;
    m_center_m_arc_end = std::atan2(tooth_tip_m.y() - m_center_m.y(), tooth_tip_m.x() - m_center_m.x());
    m_center_m_arc_end = m_center_m_arc_end < 0 ? m_center_m_arc_end + CH_C_2PI : m_center_m_arc_end;
    if (m_center_m_arc_start > m_center_m_arc_end) {
        double temp = m_center_m_arc_start;
        m_center_m_arc_start = m_center_m_arc_end;
        m_center_m_arc_end = temp;
    }

    // Express the tread body location and orientation in global frame.
    ChVector<> loc = chassis->TransformPointLocalToParent(location);
    ChQuaternion<> rot = chassis->GetRot() * rotation;
    ChVector<> xdir = rot.GetXaxis();
    ChVector<> ydir = rot.GetYaxis();
    ChVector<> zdir = rot.GetZaxis();

    // Create the tread body
    m_shoe = std::shared_ptr<ChBody>(chassis->GetSystem()->NewBody());
    m_shoe->SetNameString(m_name + "_tread");
    m_shoe->SetPos(loc);
    m_shoe->SetRot(rot);
    m_shoe->SetMass(GetTreadMass());
    m_shoe->SetInertiaXX(GetTreadInertia());
    chassis->GetSystem()->AddBody(m_shoe);

    // Add contact geometry.
    m_shoe->SetCollide(true);

    switch (m_shoe->GetContactMethod()) {
        case ChMaterialSurface::NSC:
            m_shoe->GetMaterialSurfaceNSC()->SetFriction(m_friction);
            m_shoe->GetMaterialSurfaceNSC()->SetRestitution(m_restitution);
            break;
        case ChMaterialSurface::SMC:
            m_shoe->GetMaterialSurfaceSMC()->SetFriction(m_friction);
            m_shoe->GetMaterialSurfaceSMC()->SetRestitution(m_restitution);
            m_shoe->GetMaterialSurfaceSMC()->SetYoungModulus(m_young_modulus);
            m_shoe->GetMaterialSurfaceSMC()->SetPoissonRatio(m_poisson_ratio);
            m_shoe->GetMaterialSurfaceSMC()->SetKn(m_kn);
            m_shoe->GetMaterialSurfaceSMC()->SetGn(m_gn);
            m_shoe->GetMaterialSurfaceSMC()->SetKt(m_kt);
            m_shoe->GetMaterialSurfaceSMC()->SetGt(m_gt);
            break;
    }

    AddShoeContact();

    // Create the required number of web segment bodies
    ChVector<> seg_loc = loc + (0.5 * GetToothBaseLength()) * xdir - (0.5 * GetBeltWidth()) * ydir;
    //for (int is = 0; is < GetNumWebSegments(); is++) {
    //    m_web_segments.push_back(std::shared_ptr<ChBody>(chassis->GetSystem()->NewBody()));
    //    m_web_segments[is]->SetNameString(m_name + "_web_" + std::to_string(is));
    //    m_web_segments[is]->SetPos(seg_loc + ((2 * is + 1) * m_seg_length / 2) * xdir);
    //    m_web_segments[is]->SetRot(rot);
    //    m_web_segments[is]->SetMass(m_seg_mass);
    //    m_web_segments[is]->SetInertiaXX(m_seg_inertia);
    //    chassis->GetSystem()->AddBody(m_web_segments[is]);

    //    // Add contact geometry.
    //    m_web_segments[is]->SetCollide(true);

    //    switch (m_web_segments[is]->GetContactMethod()) {
    //        case ChMaterialSurface::NSC:
    //            m_web_segments[is]->GetMaterialSurfaceNSC()->SetFriction(m_friction);
    //            m_web_segments[is]->GetMaterialSurfaceNSC()->SetRestitution(m_restitution);
    //            break;
    //        case ChMaterialSurface::SMC:
    //            m_web_segments[is]->GetMaterialSurfaceSMC()->SetFriction(m_friction);
    //            m_web_segments[is]->GetMaterialSurfaceSMC()->SetRestitution(m_restitution);
    //            m_web_segments[is]->GetMaterialSurfaceSMC()->SetYoungModulus(m_young_modulus);
    //            m_web_segments[is]->GetMaterialSurfaceSMC()->SetPoissonRatio(m_poisson_ratio);
    //            m_web_segments[is]->GetMaterialSurfaceSMC()->SetKn(m_kn);
    //            m_web_segments[is]->GetMaterialSurfaceSMC()->SetGn(m_gn);
    //            m_web_segments[is]->GetMaterialSurfaceSMC()->SetKt(m_kt);
    //            m_web_segments[is]->GetMaterialSurfaceSMC()->SetGt(m_gt);
    //            break;
    //    }

    //    AddWebContact(m_web_segments[is]);
    //}

#ifdef USE_ANCF_4
    m_web_mesh = std::make_shared<ChMesh>();

    int N_x = m_num_elements_length + 1;
    int N_y = m_num_elements_width + 1;

    double dx = GetWebLength() / m_num_elements_length;
    double dy = GetBeltWidth() / m_num_elements_width;

    double dz_steel = 0.05*25.4 / 1000.0;
    double dz_rubber = (GetWebThickness()- dz_steel)/2;
    //dz_rubber = GetWebThickness();

    // Create an orthotropic material.
    // All layers for all elements share the same material.
    double rho_rubber = 1.1e3;
    ChVector<> E_rubber(0.01e9, 0.01e9, 0.01e9);
    ChVector<> nu_rubber(0.49, 0.49, 0.49);
    //ChVector<> G(0.0003e9, 0.0003e9, 0.0003e9);
    ChVector<> G_rubber = E_rubber / (2 * (1 + .49));
    auto mat_rubber = std::make_shared<ChMaterialShellANCF>(rho_rubber, E_rubber, nu_rubber, G_rubber);

    // Create an orthotropic material.
    // All layers for all elements share the same material.
    double rho_steel = 7900.0;
    ChVector<> E_steel(210e9, 210e9, 210e9);
    ChVector<> nu_steel(0.3, 0.3, 0.3);
    //ChVector<> G(0.0003e9, 0.0003e9, 0.0003e9);
    ChVector<> G_steel = E_steel / (2 * (1 + .3));
    auto mat_steel = std::make_shared<ChMaterialShellANCF>(rho_steel, E_steel, nu_steel, G_steel);


    // Create and add the nodes
    for (int x_idx = 0; x_idx < N_x; x_idx++) {
        for (int y_idx = 0; y_idx < N_y; y_idx++) {

            // Node location
            auto node_loc = seg_loc + x_idx*dx*xdir + y_idx*dy*ydir;

            // Node direction
            auto node_dir = zdir;

            // Create the node
            auto node = std::make_shared<ChNodeFEAxyzD>(node_loc, node_dir);

            node->SetMass(0);

            // Add node to mesh
            m_web_mesh->AddNode(node);
        }
    }

    // Create the elements
    for (int x_idx = 0; x_idx < m_num_elements_length; x_idx++) {
        for (int y_idx = 0; y_idx < m_num_elements_width; y_idx++) {
            // Adjacent nodes
            int node0 = y_idx + x_idx * N_y;
            int node1 = y_idx + (x_idx + 1) * N_y;
            int node2 = (y_idx + 1) + (x_idx + 1) * N_y;
            int node3 = (y_idx + 1) + x_idx * N_y;

            // Create the element and set its nodes.
            auto element = std::make_shared<ChElementShellANCF>();
            element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzD>(m_web_mesh->GetNode(node0)),
                std::dynamic_pointer_cast<ChNodeFEAxyzD>(m_web_mesh->GetNode(node1)),
                std::dynamic_pointer_cast<ChNodeFEAxyzD>(m_web_mesh->GetNode(node2)),
                std::dynamic_pointer_cast<ChNodeFEAxyzD>(m_web_mesh->GetNode(node3)));

            // Set element dimensions
            element->SetDimensions(dx, dy);

            // Add a single layers with a fiber angle of 0 degrees.
            element->AddLayer(dz_rubber, 0 * CH_C_DEG_TO_RAD, mat_rubber);
            element->AddLayer(dz_steel, 0 * CH_C_DEG_TO_RAD, mat_steel);
            element->AddLayer(dz_rubber, 0 * CH_C_DEG_TO_RAD, mat_rubber);

            // Set other element properties
            element->SetAlphaDamp(0.05);    // Structural damping for this element
            element->SetGravityOn(false);  // turn internal gravitational force calculation off

                                            // Add element to mesh
            m_web_mesh->AddElement(element);
        }
    }

    // Add the mesh to the system
    m_shoe->GetSystem()->Add(m_web_mesh);
#endif

#ifdef USE_ANCF_8
    m_web_mesh = std::make_shared<ChMesh>();

    int N_x_edge = 2 * m_num_elements_length + 1;
    int N_y_edge = 2 * m_num_elements_width + 1;
    int N_x_mid = m_num_elements_length + 1;
    int N_y_mid = m_num_elements_width + 1;

    double dx = GetWebLength() / (2*m_num_elements_length);
    double dy = GetBeltWidth() / (2*m_num_elements_width);

    double dz_steel = 0.05*25.4 / 1000.0;
    double dz_rubber = (GetWebThickness() - dz_steel) / 2;
    //dz_rubber = GetWebThickness();

    // Create an orthotropic material.
    // All layers for all elements share the same material.
    double rho_rubber = 1.1e3;
    ChVector<> E_rubber(0.01e9, 0.01e9, 0.01e9);
    ChVector<> nu_rubber(0.49, 0.49, 0.49);
    //ChVector<> G(0.0003e9, 0.0003e9, 0.0003e9);
    ChVector<> G_rubber = E_rubber / (2 * (1 + .49));
    auto mat_rubber = std::make_shared<ChMaterialShellANCF_8>(rho_rubber, E_rubber, nu_rubber, G_rubber);

    // Create an orthotropic material.
    // All layers for all elements share the same material.
    double rho_steel = 7900.0;
    ChVector<> E_steel(210e9, 210e9, 210e9);
    ChVector<> nu_steel(0.3, 0.3, 0.3);
    //ChVector<> G(0.0003e9, 0.0003e9, 0.0003e9);
    ChVector<> G_steel = E_steel / (2 * (1 + .3));
    auto mat_steel = std::make_shared<ChMaterialShellANCF_8>(rho_steel, E_steel, nu_steel, G_steel);


    // Create and add the nodes
    for (int x_idx = 0; x_idx < N_x_edge; x_idx++) {
        for (int y_idx = 0; y_idx < N_y_edge; y_idx++) {
            if ((x_idx % 2 == 1) && (y_idx % 2 == 1))
                continue;

            // Node location
            auto node_loc = seg_loc + x_idx*dx*xdir + y_idx*dy*ydir;

            // Node direction
            auto node_dir = zdir;

            // Node direction derivative
            auto node_curv = ChVector<>(0.0, 0.0, 0.0);

            // Create the node
            auto node = std::make_shared<ChNodeFEAxyzDD>(node_loc, node_dir, node_curv);

            node->SetMass(0);

            // Add node to mesh
            m_web_mesh->AddNode(node);
        }
    }

    // Create the elements
    for (int x_idx = 0; x_idx < m_num_elements_length; x_idx++) {
        for (int y_idx = 0; y_idx < m_num_elements_width; y_idx++) {
            // Adjacent nodes
            /// The node numbering is in ccw fashion as in the following scheme:
            ///         v
            ///         ^
            /// D o-----G-----o C
            ///   |     |     |
            /// --H-----+-----F-> u
            ///   |     |     |
            /// A o-----E-----o B

            int node0 = 2 * y_idx + x_idx * (N_y_edge + N_y_mid);
            int node1 = 2 * y_idx + (x_idx + 1) * (N_y_edge + N_y_mid);
            int node2 = 2 * (y_idx + 1) + (x_idx + 1) * (N_y_edge + N_y_mid);
            int node3 = 2 * (y_idx + 1) + x_idx * (N_y_edge + N_y_mid);

            int node4 = 1 * y_idx + x_idx * (N_y_edge + N_y_mid) + N_y_edge;
            int node5 = 2 * y_idx + (x_idx + 1) * (N_y_edge + N_y_mid) + 1;
            int node6 = 1 * y_idx + x_idx * (N_y_edge + N_y_mid) + N_y_edge + 1;
            int node7 = 2 * y_idx + x_idx * (N_y_edge + N_y_mid) + 1;


            // Create the element and set its nodes.
            auto element = std::make_shared<ChElementShellANCF_8>();
            element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node0)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node1)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node2)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node3)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node4)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node5)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node6)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node7)));

            // Set element dimensions
            element->SetDimensions(2*dx, 2*dy);

            // Add a single layers with a fiber angle of 0 degrees.
            element->AddLayer(dz_rubber, 0 * CH_C_DEG_TO_RAD, mat_rubber);
            element->AddLayer(dz_steel, 0 * CH_C_DEG_TO_RAD, mat_steel);
            element->AddLayer(dz_rubber, 0 * CH_C_DEG_TO_RAD, mat_rubber);

            // Set other element properties
            element->SetAlphaDamp(0.05);    // Structural damping for this element
            element->SetGravityOn(false);  // turn internal gravitational force calculation off

                                           // Add element to mesh
            m_web_mesh->AddElement(element);
        }
    }

    // Add the mesh to the system
    m_shoe->GetSystem()->Add(m_web_mesh);
#endif


#if FALSE
    if (GetIndex() == 0) {
        //-------------------------------------------------------------------
        //Problem Geometry
        double web_angle = -30 * CH_C_DEG_TO_RAD;
        double length = 24 * 25.4 / 1000;
        double width = 12 * 25.4 / 1000;
        double thickness = 1 * 25.4 / 1000;

        int num_elements_length = 10;
        int num_elements_width = 10;
        int num_elements_thickness = 1;

        int N_x = num_elements_length + 1;
        int N_y = num_elements_width + 1;

        double dx = length * std::cos(web_angle) / num_elements_length;
        double dy = width / num_elements_width;
        double dz = length * std::sin(web_angle) / num_elements_length;

        // Number of elements in the z direction is considered as 1
        int TotalNumElements = num_elements_length * num_elements_width;
        int TotalNumNodes = N_x * N_y;

        // Create a mesh, that is a container for groups of elements and their referenced nodes.
        auto my_mesh = std::make_shared<ChMesh>();

        // Create and add the nodes
        for (int x_idx = 0; x_idx < N_x; x_idx++) {
            for (int y_idx = 0; y_idx < N_y; y_idx++) {

                // Node location
                double loc_x = x_idx * dx;
                double loc_y = y_idx * dy;
                double loc_z = x_idx * dz;

                // Node direction
                double dir_x = std::cos(web_angle + CH_C_PI_2);
                double dir_y = 0;
                double dir_z = std::sin(web_angle + CH_C_PI_2);

                // Create the node
                auto node = std::make_shared<ChNodeFEAxyzD>(ChVector<>(loc_x, loc_y, loc_z), ChVector<>(dir_x, dir_y, dir_z));

                node->SetMass(0);

                // Add node to mesh
                my_mesh->AddNode(node);
            }
        }

        // Create an orthotropic material.
        // All layers for all elements share the same material.
        double rho = 1.1e3;
        ChVector<> E(0.1e9, 0.1e9, 0.1e9);
        ChVector<> nu(0.49, 0.49, 0.49);
        //ChVector<> G(0.0003e9, 0.0003e9, 0.0003e9);
        ChVector<> G = E / (2 * (1 + .49));
        auto mat = std::make_shared<ChMaterialShellANCF>(rho, E, nu, G);


        // Create the elements
        for (int x_idx = 0; x_idx < num_elements_length; x_idx++) {
            for (int y_idx = 0; y_idx < num_elements_width; y_idx++) {
                // Adjacent nodes
                int node0 = y_idx + x_idx * N_y;
                int node1 = y_idx + (x_idx + 1) * N_y;
                int node2 = (y_idx + 1) + (x_idx + 1) * N_y;
                int node3 = (y_idx + 1) + x_idx * N_y;

                // Create the element and set its nodes.
                auto element = std::make_shared<ChElementShellANCF>();
                element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node0)),
                    std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node1)),
                    std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node2)),
                    std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node3)));

                // Set element dimensions
                element->SetDimensions(std::sqrt(dx*dx + dz*dz), dy);

                // Add a single layers with a fiber angle of 0 degrees.
                element->AddLayer(dz, 0 * CH_C_DEG_TO_RAD, mat);

                // Set other element properties
                element->SetAlphaDamp(0.05);    // Structural damping for this element
                element->SetGravityOn(false);  // turn internal gravitational force calculation off

                                               // Add element to mesh
                my_mesh->AddElement(element);
            }
        }

        // Add the mesh to the system
        chassis->GetSystem()->Add(my_mesh);

        // -------------------------------------
        // Options for visualization in irrlicht
        // -------------------------------------

        auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
        mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
        mvisualizemesh->SetColorscaleMinMax(0.0, 5.50);
        mvisualizemesh->SetShrinkElements(true, 0.85);
        mvisualizemesh->SetSmoothFaces(true);
        my_mesh->AddAsset(mvisualizemesh);

        auto mvisualizemeshref = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
        mvisualizemeshref->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
        mvisualizemeshref->SetWireframe(true);
        mvisualizemeshref->SetDrawInUndeformedReference(true);
        my_mesh->AddAsset(mvisualizemeshref);

        auto mvisualizemeshC = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
        mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);
        mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
        mvisualizemeshC->SetSymbolsThickness(0.004);
        my_mesh->AddAsset(mvisualizemeshC);

        auto mvisualizemeshD = std::make_shared<ChVisualizationFEAmesh>(*(my_mesh.get()));
        // mvisualizemeshD->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_VECT_SPEED);
        mvisualizemeshD->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_ELEM_TENS_STRAIN);
        mvisualizemeshD->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
        mvisualizemeshD->SetSymbolsScale(1);
        mvisualizemeshD->SetColorscaleMinMax(-0.5, 5);
        mvisualizemeshD->SetZbufferHide(false);
        my_mesh->AddAsset(mvisualizemeshD);
        //-------------------------------------------------------------------
    }

#endif


}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChTrackShoeRigidANCFCB::Initialize(std::shared_ptr<ChBodyAuxRef> chassis,
                                    const std::vector<ChCoordsys<>>& component_pos) {
    // Check the number of provided locations and orientations.
    assert(component_pos.size() == GetNumWebSegments() + 1);

    // Initialize at origin.
    Initialize(chassis, VNULL, QUNIT);

    // Overwrite absolute body locations and orientations.
    m_shoe->SetPos(chassis->TransformPointLocalToParent(component_pos[0].pos));
    m_shoe->SetRot(chassis->GetRot() * component_pos[0].rot);

    //for (int is = 0; is < GetNumWebSegments(); is++) {
    //    m_web_segments[is]->SetPos(chassis->TransformPointLocalToParent(component_pos[is + 1].pos));
    //    m_web_segments[is]->SetRot(chassis->GetRot() * component_pos[is + 1].rot);
    //}

    // Overwrite absolute node locations and orientations.

    auto rot = chassis->GetRot() * component_pos[1].rot;
    ChVector<> xdir = rot.GetXaxis();
    ChVector<> ydir = rot.GetYaxis();
    ChVector<> zdir = rot.GetZaxis();

    ChVector<> seg_loc = chassis->TransformPointLocalToParent(component_pos[1].pos) - (0.5 * GetWebLength()) * xdir;

#ifdef USE_ANCF_4
    int N_x = m_num_elements_length + 1;
    int N_y = m_num_elements_width + 1;

    double dx = GetWebLength() / m_num_elements_length;
    double dy = GetBeltWidth() / m_num_elements_width;
    double dz = GetWebThickness();


    // Move the nodes on the mesh to the correct location
    for (int x_idx = 0; x_idx < N_x; x_idx++) {
        for (int y_idx = 0; y_idx < N_y; y_idx++) {

            // Node location
            auto node_loc = seg_loc + x_idx*dx*xdir + y_idx*dy*ydir;

            // Node direction
            auto node_dir = zdir;

            // Create the node
            int node_idx = y_idx + x_idx * N_y;
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyzD>(m_web_mesh->GetNode(node_idx));

            node->SetPos(node_loc);
            node->SetD(node_dir);
        }
    }
#endif

#ifdef USE_ANCF_8

    int N_x_edge = 2 * m_num_elements_length + 1;
    int N_y_edge = 2 * m_num_elements_width + 1;
    int N_x_mid = m_num_elements_length + 1;
    int N_y_mid = m_num_elements_width + 1;

    double dx = GetWebLength() / (2 * m_num_elements_length);
    double dy = GetBeltWidth() / (2 * m_num_elements_width);

    // Move the nodes on the mesh to the correct location
    int node_idx = 0;
    for (int x_idx = 0; x_idx < m_num_elements_length; x_idx++) {
        for (int y_idx = 0; y_idx < m_num_elements_width; y_idx++) {
            if ((x_idx % 2 == 1) && (y_idx % 2 == 1))
                continue;

            // Node location
            auto node_loc = seg_loc + x_idx*dx*xdir + y_idx*dy*ydir;

            // Node direction
            auto node_dir = zdir;

            // Create the node
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node_idx));

            node->SetPos(node_loc);
            node->SetD(node_dir);
            node_idx++;
        }
    }
#endif

}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
double ChTrackShoeRigidANCFCB::GetMass() const {
    return GetTreadMass() + GetWebMass();
}

double ChTrackShoeRigidANCFCB::GetPitch() const {
    return GetToothBaseLength() + GetWebLength();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChTrackShoeRigidANCFCB::AddShoeContact() {
    m_shoe->GetCollisionModel()->ClearModel();

    m_shoe->GetCollisionModel()->SetFamily(TrackedCollisionFamily::SHOES);
    m_shoe->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(TrackedCollisionFamily::SHOES);

    // Guide pin
    ChVector<> g_hdims = GetGuideBoxDimensions() / 2;
    ChVector<> g_loc(GetGuideBoxOffsetX(), 0, GetWebThickness() / 2 + g_hdims.z());
    m_shoe->GetCollisionModel()->AddBox(g_hdims.x(), g_hdims.y(), g_hdims.z(), g_loc);

    // Main box
    ChVector<> b_hdims(GetToothBaseLength() / 2, GetBeltWidth() / 2, GetWebThickness() / 2);
    ChVector<> b_loc(0, 0, 0);
    m_shoe->GetCollisionModel()->AddBox(b_hdims.x(), b_hdims.y(), b_hdims.z(), b_loc);

    // Tread box
    ChVector<> t_hdims(GetTreadLength() / 2, GetBeltWidth() / 2, GetTreadThickness() / 2);
    ChVector<> t_loc(0, 0, (-GetWebThickness() - GetTreadThickness()) / 2);
    m_shoe->GetCollisionModel()->AddBox(t_hdims.x(), t_hdims.y(), t_hdims.z(), t_loc);

    m_shoe->GetCollisionModel()->BuildModel();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChTrackShoeRigidANCFCB::AddWebContact(std::shared_ptr<ChBody> segment) {
    segment->GetCollisionModel()->ClearModel();

    segment->GetCollisionModel()->SetFamily(TrackedCollisionFamily::SHOES);
    segment->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(TrackedCollisionFamily::SHOES);

    segment->GetCollisionModel()->AddBox(m_seg_length / 2, GetBeltWidth() / 2, GetWebThickness() / 2);

    segment->GetCollisionModel()->BuildModel();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChTrackShoeRigidANCFCB::AddVisualizationAssets(VisualizationType vis) {
    if (vis == VisualizationType::NONE)
        return;

    AddShoeVisualization();
    for (auto segment : m_web_segments)
        AddWebVisualization(segment);
}

void ChTrackShoeRigidANCFCB::RemoveVisualizationAssets() {
    m_shoe->GetAssets().clear();
    for (auto segment : m_web_segments) {
        segment->GetAssets().clear();
    }
}

ChColor GetColor_RigidANCF(size_t index) {
    if (index == 0)
        return ChColor(0.7f, 0.4f, 0.4f);
    else if (index % 2 == 0)
        return ChColor(0.4f, 0.7f, 0.4f);
    else
        return ChColor(0.4f, 0.4f, 0.7f);
}

void ChTrackShoeRigidANCFCB::AddShoeVisualization() {
    m_shoe->AddAsset(std::make_shared<ChColorAsset>(GetColor_RigidANCF(m_index)));

    // Guide pin
    ChVector<> g_hdims = GetGuideBoxDimensions() / 2;
    ChVector<> g_loc(GetGuideBoxOffsetX(), 0, GetWebThickness() / 2 + g_hdims.z());
    auto box_pin = std::make_shared<ChBoxShape>();
    box_pin->GetBoxGeometry().Size = g_hdims;
    box_pin->GetBoxGeometry().Pos = g_loc;
    m_shoe->AddAsset(box_pin);

    // Main box
    ChVector<> b_hdims(GetToothBaseLength() / 2, GetBeltWidth() / 2, GetWebThickness() / 2);
    ChVector<> b_loc(0, 0, 0);
    auto box_main = std::make_shared<ChBoxShape>();
    box_main->GetBoxGeometry().Size = b_hdims;
    box_main->GetBoxGeometry().Pos = b_loc;
    m_shoe->AddAsset(box_main);

    // Tread box
    ChVector<> t_hdims(GetTreadLength() / 2, GetBeltWidth() / 2, GetTreadThickness() / 2);
    ChVector<> t_loc(0, 0, (-GetWebThickness() - GetTreadThickness()) / 2);
    auto box_tread = std::make_shared<ChBoxShape>();
    box_tread->GetBoxGeometry().Size = t_hdims;
    box_tread->GetBoxGeometry().Pos = t_loc;
    m_shoe->AddAsset(box_tread);

    // Connection to first web segment
    double radius = GetWebThickness() / 4;
    auto cyl = std::make_shared<ChCylinderShape>();
    cyl->GetCylinderGeometry().rad = radius;
    cyl->GetCylinderGeometry().p1 = ChVector<>(GetToothBaseLength() / 2, -GetBeltWidth() / 2 - 2 * radius, 0);
    cyl->GetCylinderGeometry().p2 = ChVector<>(GetToothBaseLength() / 2, +GetBeltWidth() / 2 + 2 * radius, 0);
    m_shoe->AddAsset(cyl);

    // Create tooth meshes
    m_shoe->AddAsset(ToothMesh(GetBeltWidth() / 2 - GetToothWidth() / 2));
    m_shoe->AddAsset(ToothMesh(-GetBeltWidth() / 2 + GetToothWidth() / 2));
}

void ChTrackShoeRigidANCFCB::AddWebVisualization(std::shared_ptr<ChBody> segment) {
    segment->AddAsset(std::make_shared<ChColorAsset>(GetColor_RigidANCF(m_index)));

    //auto box = std::make_shared<ChBoxShape>();
    //box->GetBoxGeometry().SetLengths(ChVector<>(m_seg_length, GetBeltWidth(), GetWebThickness()));
    //segment->AddAsset(box);

    auto cyl = std::make_shared<ChCylinderShape>();
    double radius = GetWebThickness() / 4;
    cyl->GetCylinderGeometry().rad = radius;
    cyl->GetCylinderGeometry().p1 = ChVector<>(m_seg_length / 2, -GetBeltWidth() / 2 - 2 * radius, 0);
    cyl->GetCylinderGeometry().p2 = ChVector<>(m_seg_length / 2, +GetBeltWidth() / 2 + 2 * radius, 0);
    segment->AddAsset(cyl);

    auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*(m_web_mesh.get()));
    mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
    mvisualizemesh->SetColorscaleMinMax(0.0, 5.50);
    mvisualizemesh->SetShrinkElements(true, 0.85);
    mvisualizemesh->SetSmoothFaces(true);
    m_web_mesh->AddAsset(mvisualizemesh);

    auto mvisualizemeshref = std::make_shared<ChVisualizationFEAmesh>(*(m_web_mesh.get()));
    mvisualizemeshref->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizemeshref->SetWireframe(true);
    mvisualizemeshref->SetDrawInUndeformedReference(true);
    m_web_mesh->AddAsset(mvisualizemeshref);

    auto mvisualizemeshC = std::make_shared<ChVisualizationFEAmesh>(*(m_web_mesh.get()));
    mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);
    mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizemeshC->SetSymbolsThickness(0.004);
    m_web_mesh->AddAsset(mvisualizemeshC);

    auto mvisualizemeshD = std::make_shared<ChVisualizationFEAmesh>(*(m_web_mesh.get()));
    // mvisualizemeshD->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_VECT_SPEED);
    mvisualizemeshD->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_ELEM_TENS_STRAIN);
    mvisualizemeshD->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizemeshD->SetSymbolsScale(1);
    mvisualizemeshD->SetColorscaleMinMax(-0.5, 5);
    mvisualizemeshD->SetZbufferHide(false);
    m_web_mesh->AddAsset(mvisualizemeshD);

}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChTrackShoeRigidANCFCB::Connect(std::shared_ptr<ChTrackShoe> next) {
    ChSystem* system = m_shoe->GetSystem();
    ChVector<> loc_cur_shoe = m_shoe->TransformPointLocalToParent(ChVector<>(GetToothBaseLength() / 2, 0, 0));
    ChQuaternion<> rot_cur_shoe = m_shoe->GetRot();
    ChVector<> loc_next_shoe = next->GetShoeBody()->TransformPointLocalToParent(ChVector<>(-GetToothBaseLength() / 2, 0, 0));;
    ChQuaternion<> rot_next_shoe = next->GetShoeBody()->GetRot();

#ifdef USE_ANCF_4
    int N_x = m_num_elements_length + 1;
    int N_y = m_num_elements_width + 1;

    // Change the gradient on the web boundary nodes that will connect to the current shoe body
    // and then connect those web nodes to the show tread body
    for (int y_idx = 0; y_idx < m_num_elements_width; y_idx++) {
        int node_idx = y_idx + 0 * N_y;
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzD>(m_web_mesh->GetNode(node_idx));

        node->SetD(rot_cur_shoe.GetZaxis());

        auto constraintxyz = std::make_shared<ChLinkPointFrame>();
        constraintxyz->Initialize(node, m_shoe);
        system->Add(constraintxyz);

        auto constraintD = std::make_shared<ChLinkDirFrame>();
        constraintD->Initialize(node, m_shoe);
        system->Add(constraintD);
    }

    // Change the gradient on the boundary nodes that will connect to the second fixed body
    // and then connect those nodes to the body
    for (int y_idx = 0; y_idx < m_num_elements_width; y_idx++) {
        int node_idx = y_idx + m_num_elements_length * N_y;
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzD>(m_web_mesh->GetNode(node_idx));

        node->SetD(rot_next_shoe.GetZaxis());

        auto constraintxyz = std::make_shared<ChLinkPointFrame>();
        constraintxyz->Initialize(node, next->GetShoeBody());
        system->Add(constraintxyz);

        auto constraintD = std::make_shared<ChLinkDirFrame>();
        constraintD->Initialize(node, next->GetShoeBody());
        system->Add(constraintD);
    }
#endif

#ifdef USE_ANCF_8
    int N_x_edge = 2 * m_num_elements_length + 1;
    int N_y_edge = 2 * m_num_elements_width + 1;
    int N_x_mid = m_num_elements_length + 1;
    int N_y_mid = m_num_elements_width + 1;

    double dx = GetWebLength() / (2 * m_num_elements_length);
    double dy = GetBeltWidth() / (2 * m_num_elements_width);

    // Change the gradient on the web boundary nodes that will connect to the current shoe body
    // and then connect those web nodes to the show tread body
    for (int y_idx = 0; y_idx < N_y_edge; y_idx++) {
        int node_idx = y_idx;
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node_idx));

        node->SetD(rot_cur_shoe.GetZaxis());

        auto constraintxyz = std::make_shared<ChLinkPointFrame>();
        constraintxyz->Initialize(node, m_shoe);
        system->Add(constraintxyz);

        auto constraintD = std::make_shared<ChLinkDirFrame>();
        constraintD->Initialize(node, m_shoe);
        system->Add(constraintD);
    }

    // Change the gradient on the boundary nodes that will connect to the second fixed body
    // and then connect those nodes to the body
    for (int y_idx = 0; y_idx < N_y_edge; y_idx++) {
        int node_idx = y_idx + m_num_elements_length * (N_y_edge + N_y_mid);
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzDD>(m_web_mesh->GetNode(node_idx));

        node->SetD(rot_next_shoe.GetZaxis());

        auto constraintxyz = std::make_shared<ChLinkPointFrame>();
        constraintxyz->Initialize(node, next->GetShoeBody());
        system->Add(constraintxyz);

        auto constraintD = std::make_shared<ChLinkDirFrame>();
        constraintD->Initialize(node, next->GetShoeBody());
        system->Add(constraintD);
    }
#endif

}


// -----------------------------------------------------------------------------
// Utilities for creating tooth mesh
// -----------------------------------------------------------------------------
size_t ChTrackShoeRigidANCFCB::ProfilePoints(std::vector<ChVector2<>>& points, std::vector<ChVector2<>>& normals) {
    int np = 4;
    double step = 1.0 / (np - 1);

    // Start from point on left tooth base (Am).
    ChVector2<> Am(-GetToothBaseLength() / 2, GetWebThickness() / 2);
    ChVector2<> Bm(-GetToothTipLength() / 2, GetToothHeight() + GetWebThickness() / 2);
    for (int i = 0; i < np; i++) {
        ChVector2<> pc(Am.x() + (i * step) * (Bm.x() - Am.x()), Am.y() + (i * step) * (Bm.y() - Am.y()));
        ChVector2<> nrm = (pc - m_center_m).GetNormalized();
        ChVector2<> pt = m_center_m + nrm * GetToothArcRadius();
        points.push_back(pt);
        normals.push_back(nrm);
    }

    // Mid-point on tooth tip.
    points.push_back(ChVector2<>(0, GetToothHeight() + GetWebThickness() / 2));
    normals.push_back(ChVector2<>(0, 1));

    // Continue from point on right tooth tip (Bp).
    ChVector2<> Ap(GetToothBaseLength() / 2, GetWebThickness() / 2);
    ChVector2<> Bp(GetToothTipLength() / 2, GetToothHeight() + GetWebThickness() / 2);
    for (int i = 0; i < np; i++) {
        ChVector2<> pc(Bp.x() + (i * step) * (Ap.x() - Bp.x()), Bp.y() + (i * step) * (Ap.y() - Bp.y()));
        ChVector2<> nrm = (pc - m_center_p).GetNormalized();
        ChVector2<> pt = m_center_p + nrm * GetToothArcRadius();
        points.push_back(pt);
        normals.push_back(nrm);
    }

    ////std::cout << std::endl << std::endl;
    ////for (auto p : points)
    ////    std::cout << p.x() << "  " << p.y() << std::endl;

    return points.size();
}

std::shared_ptr<ChTriangleMeshShape> ChTrackShoeRigidANCFCB::ToothMesh(double y) {
    // Obtain profile points.
    std::vector<ChVector2<>> points2;
    std::vector<ChVector2<>> normals2;
    size_t np = ProfilePoints(points2, normals2);

    // Create the triangular mesh.
    geometry::ChTriangleMeshConnected trimesh;
    std::vector<ChVector<>>& vertices = trimesh.getCoordsVertices();
    std::vector<ChVector<>>& normals = trimesh.getCoordsNormals();
    std::vector<ChVector<int>>& idx_vertices = trimesh.getIndicesVertexes();
    std::vector<ChVector<int>>& idx_normals = trimesh.getIndicesNormals();

    // Number of vertices:
    //   - 1 for the middle of the tooth base on +y side
    //   - np for the +y tooth side
    //   - 1 for the middle of the tooth base on -y side
    //   - np for the -y tooth side
    size_t num_vertices = 2 * (np + 1);
    vertices.resize(num_vertices);

    // Number of normals:
    //   - 1 for the +y face
    //   - 1 for the -y face
    //   - np for the tooth surface
    size_t num_normals = 2 + np;
    normals.resize(num_normals);

    // Number of faces:
    //    - np-1 for the +y side
    //    - np-1 for the -y side
    //    - 2 * (np-1) for the tooth surface
    size_t num_faces = 4 * (np - 1);
    idx_vertices.resize(num_faces);
    idx_normals.resize(num_faces);

    // Load vertices.
    double yp = y + GetToothWidth() / 2;
    double ym = y - GetToothWidth() / 2;

    size_t iv = 0;
    vertices[iv++] = ChVector<>(0, yp, GetWebThickness() / 2);
    for (size_t i = 0; i < np; i++)
        vertices[iv++] = ChVector<>(points2[i].x(), yp, points2[i].y());
    vertices[iv++] = ChVector<>(0, ym, GetWebThickness() / 2);
    for (size_t i = 0; i < np; i++)
        vertices[iv++] = ChVector<>(points2[i].x(), ym, points2[i].y());

    // Load normals.
    size_t in = 0;
    normals[in++] = ChVector<>(0, +1, 0);
    normals[in++] = ChVector<>(0, -1, 0);
    for (size_t i = 0; i < np; i++)
        normals[in++] = ChVector<>(normals2[i].x(), 0, normals2[i].y());

    // Load triangles on +y side.
    size_t it = 0;
    for (size_t i = 0; i < np - 1; i++) {
        idx_vertices[it] = ChVector<int>(0, i + 1, i + 2);
        idx_normals[it] = ChVector<int>(0, 0, 0);
        it++;
    }

    // Load triangles on -y side.
    for (size_t i = 0; i < np - 1; i++) {
        idx_vertices[it] = ChVector<int>(0, i + 1, i + 2) + (np + 1);
        idx_normals[it] = ChVector<int>(1, 1, 1);
        it++;
    }

    // Load triangles on tooth surface.
    for (size_t i = 0; i < np - 1; i++) {
        idx_vertices[it] = ChVector<int>(i + 1, i + 1 + (np + 1), i + 2 + (np + 1));
        idx_normals[it] = ChVector<int>(i + 2, i + 2, i + 3);
        it++;
        idx_vertices[it] = ChVector<int>(i + 1, i + 2 + (np + 1), i + 2);
        idx_normals[it] = ChVector<int>(i + 2, i + 3, i + 3);
        it++;
    }

    auto trimesh_shape = std::make_shared<ChTriangleMeshShape>();
    trimesh_shape->SetMesh(trimesh);

    return trimesh_shape;
}

}  // end namespace vehicle
}  // end namespace chrono
