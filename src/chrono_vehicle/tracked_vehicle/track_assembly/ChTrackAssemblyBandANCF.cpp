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
// Base class for a continuous band track assembly using an ANCFshell-based web
// (template definition).
//
// The reference frame for a vehicle follows the ISO standard: Z-axis up, X-axis
// pointing forward, and Y-axis towards the left of the vehicle.
//
// =============================================================================

#include <cmath>

#include "chrono/core/ChLog.h"
#include "chrono_vehicle/tracked_vehicle/track_assembly/ChTrackAssemblyBandANCF.h"

using namespace chrono::fea;

namespace chrono {
namespace vehicle {

// -----------------------------------------------------------------------------
// Assemble the track shoes over the wheels.
//
// Returns true if the track shoes were initialized in a counter clockwise
// direction and false otherwise.
//
// This procedure is performed in the chassis reference frame, taking into
// account the convention that the chassis reference frame has the x-axis
// pointing to the front of the vehicle and the z-axis pointing up.
// It is also assumed that the sprocket, idler, and road wheels lie in the
// same vertical plane (in the chassis reference frame). The assembly is done
// in the (x-z) plane.
//
// TODO: NEEDS fixes for clock-wise wrapping (idler in front of sprocket)
//
// -----------------------------------------------------------------------------
bool ChTrackAssemblyBandANCF::Assemble(std::shared_ptr<ChBodyAuxRef> chassis) {
    // Number of track shoes
    int num_shoes = static_cast<int>(m_shoes.size());

    // Set up web connection lengths
    ChVectorDynamic<> ShoeConnectionLengths(1 + m_shoes[0]->GetNumWebSegments());
    ShoeConnectionLengths(0) = m_shoes[0]->GetToothBaseLength();

    double seg_length = m_shoes[0]->GetWebLength() / m_shoes[0]->GetNumWebSegments();
    for (int is = 0; is < m_shoes[0]->GetNumWebSegments(); is++) {
        ShoeConnectionLengths(1 + is) = seg_length;
    }

    // Calculate assembly points
    ChMatrixDynamic<> ShoePoints;
    bool ccw = FindAssemblyPoints(chassis, num_shoes, ShoeConnectionLengths, ShoePoints);

    // Create and add the mesh container for the track shoe webs to the system
    m_track_mesh = std::make_shared<ChMesh>();
    chassis->GetSystem()->Add(m_track_mesh);

    // Now create all of the track shoes at the located points
    auto num_shoe_elements = ShoeConnectionLengths.GetRows();
    for (int s = 0; s < num_shoes; s++) {
        std::vector<ChCoordsys<>> shoe_components_coordsys;
        for (int i = 0; i < num_shoe_elements; i++) {
            ChVector<> loc(
                (ShoePoints(i + s * num_shoe_elements, 0) + ShoePoints(i + 1 + s * num_shoe_elements, 0)) / 2,
                m_sprocket_offset,
                (ShoePoints(i + s * num_shoe_elements, 1) + ShoePoints(i + 1 + s * num_shoe_elements, 1)) / 2);

            double ang =
                std::atan2(ShoePoints(i + 1 + s * num_shoe_elements, 1) - ShoePoints(i + s * num_shoe_elements, 1),
                           ShoePoints(i + 1 + s * num_shoe_elements, 0) - ShoePoints(i + s * num_shoe_elements, 0));
            ChQuaternion<> rot = Q_from_AngY(-ang);  // Negative of the angle in 3D

            shoe_components_coordsys.push_back(ChCoordsys<>(loc, rot));
        }

        // Set index within the track assembly
        m_shoes[s]->SetIndex(s);
        // Pass the track mesh container to the shoe so that it adds to it
        m_shoes[s]->SetMesh(m_track_mesh);
        // Initialize the track shoe system
        m_shoes[s]->Initialize(chassis, shoe_components_coordsys);
    }

    // Add contact for the mesh
    //// TODO

    GetLog() << "Track assembly done.  Number of track shoes: " << ShoePoints.GetRows() / 2 << "\n";

    return ccw;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

void ChTrackAssemblyBandANCF::AddVisualizationAssets(VisualizationType vis) {
    if (vis == VisualizationType::NONE)
        return;

    auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*(m_track_mesh.get()));
    mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
    mvisualizemesh->SetColorscaleMinMax(0.0, 5.50);
    mvisualizemesh->SetShrinkElements(true, 0.85);
    mvisualizemesh->SetSmoothFaces(true);
    m_track_mesh->AddAsset(mvisualizemesh);

    auto mvisualizemeshref = std::make_shared<ChVisualizationFEAmesh>(*(m_track_mesh.get()));
    mvisualizemeshref->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizemeshref->SetWireframe(true);
    mvisualizemeshref->SetDrawInUndeformedReference(true);
    m_track_mesh->AddAsset(mvisualizemeshref);

    auto mvisualizemeshC = std::make_shared<ChVisualizationFEAmesh>(*(m_track_mesh.get()));
    mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);
    mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizemeshC->SetSymbolsThickness(0.004);
    m_track_mesh->AddAsset(mvisualizemeshC);
}

void ChTrackAssemblyBandANCF::RemoveVisualizationAssets() {
    m_track_mesh->GetAssets().clear();
}

}  // end namespace vehicle
}  // end namespace chrono
