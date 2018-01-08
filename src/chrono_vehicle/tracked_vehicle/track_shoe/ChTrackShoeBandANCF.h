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
// Base class for a continuous band track shoe using an ANCFshell-based web
// (template definition).
//
// =============================================================================

#ifndef CH_TRACK_SHOE_BAND_ANCF_H
#define CH_TRACK_SHOE_BAND_ANCF_H

#include "chrono_vehicle/tracked_vehicle/track_shoe/ChTrackShoeBand.h"

#include "chrono_fea/ChContactSurfaceMesh.h"
#include "chrono_fea/ChContactSurfaceNodeCloud.h"
#include "chrono_fea/ChElementShellANCF.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChNodeFEAbase.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"

namespace chrono {
namespace vehicle {

/// @addtogroup vehicle_tracked_shoe
/// @{

/// Base class for a continuous band track shoe using an ANCFshell-based web.
/// (template definition)
class CH_VEHICLE_API ChTrackShoeBandANCF : public ChTrackShoeBand {
  public:
    ChTrackShoeBandANCF(const std::string& name  ///< [in] name of the subsystem
    );

    virtual ~ChTrackShoeBandANCF() {}

    /// Initialize this track shoe subsystem.
    /// The track shoe is created within the specified system and initialized
    /// at the specified location and orientation (expressed in the global frame).
    /// This version initializes the bodies of a CB rigid-link track shoe such that
    /// the center of the track shoe subsystem is at the specified location and all
    /// bodies have the specified orientation.
    virtual void Initialize(std::shared_ptr<ChBodyAuxRef> chassis,  ///< [in] handle to the chassis body
                            const ChVector<>& location,             ///< [in] location relative to the chassis frame
                            const ChQuaternion<>& rotation          ///< [in] orientation relative to the chassis frame
                            ) override;

    /// Initialize this track shoe system.
    /// This version specifies the locations and orientations of the tread body and of
    /// the web link bodies (relative to the chassis frame).
    void Initialize(std::shared_ptr<ChBodyAuxRef> chassis,          ///< [in] handle to chassis body
                    const std::vector<ChCoordsys<>>& component_pos  ///< [in] location & orientation of the shoe bodies
    );

    /// Connect this track shoe to the specified neighbor.
    /// This function must be called only after both track shoes have been initialized.
    virtual void Connect(std::shared_ptr<ChTrackShoe> next  ///< [in] handle to the neighbor track shoe
                         ) override;

    /// Add visualization assets for the track shoe subsystem.
    virtual void AddVisualizationAssets(VisualizationType vis) override;

    /// Remove visualization assets for the track shoe subsystem.
    virtual void RemoveVisualizationAssets() override final;

  protected:
    /// Return the number of segments that the web section is broken up into
    virtual int GetNumWebSegments() const = 0;

    /// Add contact geometry for the web.
    /// Note that this is for contact with wheels, idler, and ground only.
    /// This contact geometry does not affect contact with the sprocket.
    virtual void AddWebContact();

    /// Set the FEA Mesh container to use (Can only be called before the shoe is initialized)
    /// Returns true if the FEA Mesh container was set
    /// Returns false if it was not set.
    virtual bool SetMesh(std::shared_ptr<fea::ChMesh> mesh);

  private:
    /// Add visualization of the web.
    void AddWebVisualization();

    std::shared_ptr<fea::ChMesh> m_web_mesh;

    int m_num_elements_length = 3;
    int m_num_elements_width = 4;

    bool m_is_initialized = false;

    unsigned int m_starting_node_index;

    friend class ChSprocketBandANCF;
    friend class SprocketBandANCFContactCB;
    friend class ChTrackAssemblyBandANCF;
};

/// Vector of handles to continuous band ANCFshell-based track shoe subsystems.
typedef std::vector<std::shared_ptr<ChTrackShoeBandANCF>> ChTrackShoeBandANCFList;

/// @} vehicle_tracked_shoe

}  // end namespace vehicle
}  // end namespace chrono

#endif
