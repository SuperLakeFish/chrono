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
// Base class for a continuous band track shoe using a bushing-based web
// (template definition).
//
// =============================================================================

#include "chrono/assets/ChBoxShape.h"
#include "chrono/assets/ChColorAsset.h"
#include "chrono/assets/ChCylinderShape.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/physics/ChGlobal.h"

#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"

#include "chrono_vehicle/ChSubsysDefs.h"
#include "chrono_vehicle/tracked_vehicle/track_shoe/ChTrackShoeBandBushing.h"

namespace chrono {
namespace vehicle {

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

ChTrackShoeBandBushing::ChTrackShoeBandBushing(const std::string& name) : ChTrackShoeBand(name) {}

void ChTrackShoeBandBushing::Initialize(std::shared_ptr<ChBodyAuxRef> chassis,
                                        const ChVector<>& location,
                                        const ChQuaternion<>& rotation) {
    // Initialize base class (create tread body)
    ChTrackShoeBand::Initialize(chassis, location, rotation);

    // Cache values calculated from template parameters.
    m_seg_length = GetWebLength() / GetNumWebSegments();
    m_seg_mass = GetWebMass() / GetNumWebSegments();
    m_seg_inertia = GetWebInertia();  //// TODO - properly distribute web inertia

    // Express the tread body location and orientation in global frame.
    ChVector<> loc = chassis->TransformPointLocalToParent(location);
    ChQuaternion<> rot = chassis->GetRot() * rotation;
    ChVector<> xdir = rot.GetXaxis();

    // Create the required number of web segment bodies
    ChVector<> seg_loc = loc + (0.5 * GetToothBaseLength()) * xdir;
    for (int is = 0; is < GetNumWebSegments(); is++) {
        m_web_segments.push_back(std::shared_ptr<ChBody>(chassis->GetSystem()->NewBody()));
        m_web_segments[is]->SetNameString(m_name + "_web_" + std::to_string(is));
        m_web_segments[is]->SetPos(seg_loc + ((2 * is + 1) * m_seg_length / 2) * xdir);
        m_web_segments[is]->SetRot(rot);
        m_web_segments[is]->SetMass(m_seg_mass);
        m_web_segments[is]->SetInertiaXX(m_seg_inertia);
        chassis->GetSystem()->AddBody(m_web_segments[is]);

        // Add contact geometry.
        m_web_segments[is]->SetCollide(true);

        switch (m_web_segments[is]->GetContactMethod()) {
            case ChMaterialSurface::NSC:
                m_web_segments[is]->GetMaterialSurfaceNSC()->SetFriction(m_friction);
                m_web_segments[is]->GetMaterialSurfaceNSC()->SetRestitution(m_restitution);
                break;
            case ChMaterialSurface::SMC:
                m_web_segments[is]->GetMaterialSurfaceSMC()->SetFriction(m_friction);
                m_web_segments[is]->GetMaterialSurfaceSMC()->SetRestitution(m_restitution);
                m_web_segments[is]->GetMaterialSurfaceSMC()->SetYoungModulus(m_young_modulus);
                m_web_segments[is]->GetMaterialSurfaceSMC()->SetPoissonRatio(m_poisson_ratio);
                m_web_segments[is]->GetMaterialSurfaceSMC()->SetKn(m_kn);
                m_web_segments[is]->GetMaterialSurfaceSMC()->SetGn(m_gn);
                m_web_segments[is]->GetMaterialSurfaceSMC()->SetKt(m_kt);
                m_web_segments[is]->GetMaterialSurfaceSMC()->SetGt(m_gt);
                break;
        }

        AddWebContact(m_web_segments[is]);
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChTrackShoeBandBushing::Initialize(std::shared_ptr<ChBodyAuxRef> chassis,
                                        const std::vector<ChCoordsys<>>& component_pos) {
    // Check the number of provided locations and orientations.
    assert(component_pos.size() == GetNumWebSegments() + 1);

    // Initialize at origin.
    Initialize(chassis, VNULL, QUNIT);

    // Overwrite absolute body locations and orientations.
    m_shoe->SetPos(chassis->TransformPointLocalToParent(component_pos[0].pos));
    m_shoe->SetRot(chassis->GetRot() * component_pos[0].rot);

    for (int is = 0; is < GetNumWebSegments(); is++) {
        m_web_segments[is]->SetPos(chassis->TransformPointLocalToParent(component_pos[is + 1].pos));
        m_web_segments[is]->SetRot(chassis->GetRot() * component_pos[is + 1].rot);
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChTrackShoeBandBushing::AddWebContact(std::shared_ptr<ChBody> segment) {
    segment->GetCollisionModel()->ClearModel();

    segment->GetCollisionModel()->SetFamily(TrackedCollisionFamily::SHOES);
    segment->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(TrackedCollisionFamily::SHOES);

    segment->GetCollisionModel()->AddBox(m_seg_length / 2, GetBeltWidth() / 2, GetWebThickness() / 2);

    segment->GetCollisionModel()->BuildModel();
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChTrackShoeBandBushing::AddVisualizationAssets(VisualizationType vis) {
    if (vis == VisualizationType::NONE)
        return;

    AddShoeVisualization();
    for (auto segment : m_web_segments)
        AddWebVisualization(segment);
}

void ChTrackShoeBandBushing::RemoveVisualizationAssets() {
    m_shoe->GetAssets().clear();
    for (auto segment : m_web_segments) {
        segment->GetAssets().clear();
    }
}

void ChTrackShoeBandBushing::AddWebVisualization(std::shared_ptr<ChBody> segment) {
    segment->AddAsset(std::make_shared<ChColorAsset>(GetColor(m_index)));

    auto box = std::make_shared<ChBoxShape>();
    box->GetBoxGeometry().SetLengths(ChVector<>(m_seg_length, GetBeltWidth(), GetWebThickness()));
    segment->AddAsset(box);

    auto cyl = std::make_shared<ChCylinderShape>();
    double radius = GetWebThickness() / 4;
    cyl->GetCylinderGeometry().rad = radius;
    cyl->GetCylinderGeometry().p1 = ChVector<>(m_seg_length / 2, -GetBeltWidth() / 2 - 2 * radius, 0);
    cyl->GetCylinderGeometry().p2 = ChVector<>(m_seg_length / 2, +GetBeltWidth() / 2 + 2 * radius, 0);
    segment->AddAsset(cyl);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChTrackShoeBandBushing::Connect(std::shared_ptr<ChTrackShoe> next) {
    ChSystem* system = m_shoe->GetSystem();
    ChVector<> loc;
    ChQuaternion<> rot;

#if TRUE  // Use busing elements to connect the belt segments, otherwise use revolute joints
    // Bushings are inherited from ChLoad, so they require a 'load container'

    auto my_loadcontainer = std::make_shared<ChLoadContainer>();
    system->Add(my_loadcontainer);

    ChMatrixNM<double, 6, 6> K_matrix;
    ChMatrixNM<double, 6, 6> R_matrix;

    // Sample Stiffness and Damping matrix values for testing purposes
    for (unsigned int ii = 0; ii < 3; ii++) {
        K_matrix(ii, ii) = 70000000.0;
        R_matrix(ii, ii) = 0.05 * K_matrix(ii, ii);
    }
    for (unsigned int ii = 3; ii < 6; ii++) {
        K_matrix(ii, ii) = 100000.0;
        R_matrix(ii, ii) = 0.05 * K_matrix(ii, ii);
    }
    K_matrix(4, 4) = 500;
    R_matrix(4, 4) = 0.05 * K_matrix(4, 4);

    // Connect tread body to the first web segment.
    loc = m_shoe->TransformPointLocalToParent(ChVector<>(GetToothBaseLength() / 2, 0, 0));
    rot = m_shoe->GetRot();
    auto my_loadbushingg0 = std::make_shared<ChLoadBodyBodyBushingGeneric>(
        m_shoe,               // body A
        m_web_segments[0],    // body B
        ChFrame<>(loc, rot),  // initial frame of bushing in abs space
        K_matrix,             // the 6x6 (translation+rotation) K matrix in local frame
        R_matrix              // the 6x6 (translation+rotation) R matrix in local frame
    );
    my_loadbushingg0->SetApplicationFrameA(ChFrame<>(ChVector<>(GetToothBaseLength() / 2, 0, 0)));
    my_loadbushingg0->SetApplicationFrameB(ChFrame<>(ChVector<>(-m_seg_length / 2, 0, 0)));
    my_loadcontainer->Add(my_loadbushingg0);

    // Connect the web segments to each other.
    for (size_t is = 0; is < GetNumWebSegments() - 1; is++) {
        loc = m_web_segments[is]->TransformPointLocalToParent(ChVector<>(m_seg_length / 2, 0, 0));
        rot = m_web_segments[is]->GetRot();
        auto my_loadbushingg = std::make_shared<ChLoadBodyBodyBushingGeneric>(
            m_web_segments[is],      // body A
            m_web_segments[is + 1],  // body B
            ChFrame<>(loc, rot),     // initial frame of bushing in abs space
            K_matrix,                // the 6x6 (translation+rotation) K matrix in local frame
            R_matrix                 // the 6x6 (translation+rotation) R matrix in local frame
        );
        my_loadbushingg->SetApplicationFrameA(ChFrame<>(ChVector<>(m_seg_length / 2, 0, 0)));
        my_loadbushingg->SetApplicationFrameB(ChFrame<>(ChVector<>(-m_seg_length / 2, 0, 0)));
        my_loadcontainer->Add(my_loadbushingg);
    }

    // Connect the last web segment to the tread body from the next track shoe.
    int is = GetNumWebSegments() - 1;
    loc = m_web_segments[is]->TransformPointLocalToParent(ChVector<>(m_seg_length / 2, 0, 0));
    rot = m_web_segments[is]->GetRot();
    auto my_loadbushingg1 = std::make_shared<ChLoadBodyBodyBushingGeneric>(
        m_web_segments[is],   // body A
        next->GetShoeBody(),  // body B
        ChFrame<>(loc, rot),  // initial frame of bushing in abs space
        K_matrix,             // the 6x6 (translation+rotation) K matrix in local frame
        R_matrix              // the 6x6 (translation+rotation) R matrix in local frame
    );
    my_loadbushingg1->SetApplicationFrameA(ChFrame<>(ChVector<>(m_seg_length / 2, 0, 0)));
    my_loadbushingg1->SetApplicationFrameB(ChFrame<>(ChVector<>(-GetToothBaseLength() / 2, 0, 0)));
    my_loadcontainer->Add(my_loadbushingg1);

#else

    if (m_index == 0) {
        // Create and initialize a point-line joint (sliding line along X)
        loc = m_shoe->TransformPointLocalToParent(ChVector<>(GetToothBaseLength() / 2, 0, 0));
        rot = m_shoe->GetRot() * Q_from_AngZ(CH_C_PI_2);

        auto pointline = std::make_shared<ChLinkLockPointLine>();
        pointline->SetNameString(m_name + "_pointline");
        pointline->Initialize(m_shoe, m_web_segments[0], ChCoordsys<>(loc, rot));
        system->AddLink(pointline);
    } else {
        // Connect tread body to the first web segment Via a Revolute Joint
        loc = m_shoe->TransformPointLocalToParent(ChVector<>(GetToothBaseLength() / 2, 0, 0));
        rot = m_shoe->GetRot() * Q_from_AngX(CH_C_PI_2);
        auto revolute0 = std::make_shared<ChLinkLockRevolute>();
        system->AddLink(revolute0);
        revolute0->SetNameString(m_name + "_revolute_0");
        revolute0->Initialize(m_shoe, m_web_segments[0], ChCoordsys<>(loc, rot));
    }

    // Connect the web segments to each other.
    for (size_t is = 0; is < GetNumWebSegments() - 1; is++) {
        loc = m_web_segments[is]->TransformPointLocalToParent(ChVector<>(m_seg_length / 2, 0, 0));
        rot = m_web_segments[is]->GetRot() * Q_from_AngX(CH_C_PI_2);
        auto revolute = std::make_shared<ChLinkLockRevolute>();
        system->AddLink(revolute);
        revolute->SetNameString(m_name + "_revolute_" + std::to_string(is + 1));
        revolute->Initialize(m_web_segments[is], m_web_segments[is + 1], ChCoordsys<>(loc, rot));
    }

    // Connect the last web segment to the tread body from the next track shoe.
    int is = GetNumWebSegments() - 1;
    loc = m_web_segments[is]->TransformPointLocalToParent(ChVector<>(m_seg_length / 2, 0, 0));
    rot = m_web_segments[is]->GetRot() * Q_from_AngX(CH_C_PI_2);
    auto revolute1 = std::make_shared<ChLinkLockRevolute>();
    system->AddLink(revolute1);
    revolute1->SetNameString(m_name + "_revolute");
    revolute1->Initialize(m_web_segments[is], next->GetShoeBody(), ChCoordsys<>(loc, rot));

#endif
}

}  // end namespace vehicle
}  // end namespace chrono
