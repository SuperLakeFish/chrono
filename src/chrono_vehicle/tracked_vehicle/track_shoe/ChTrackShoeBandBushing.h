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

#ifndef CH_TRACK_SHOE_BAND_BUSHING_H
#define CH_TRACK_SHOE_BAND_BUSHING_H

#include "chrono_vehicle/tracked_vehicle/track_shoe/ChTrackShoeBand.h"

namespace chrono {
namespace vehicle {

/// @addtogroup vehicle_tracked_shoe
/// @{

/// Base class for a continuous band track shoe using a bushing-based web.
/// (template definition)
class CH_VEHICLE_API ChTrackShoeBandBushing : public ChTrackShoeBand {
  public:
    ChTrackShoeBandBushing(const std::string& name  ///< [in] name of the subsystem
    );

    virtual ~ChTrackShoeBandBushing() {}

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
    /// Return the length of just one of the web sections (in the X direction)
    double GetGetWebSegmentLength() { return m_seg_length; };
    /// Return a pointer to the web segment body with the provided index
    std::shared_ptr<ChBody> GetGetWebSegment(size_t index) { return m_web_segments[index]; };

    /// Add contact geometry for a web segment body.
    /// Note that this is for contact with wheels, idler, and ground only.
    /// This contact geometry does not affect contact with the sprocket.
    virtual void AddWebContact(std::shared_ptr<ChBody> segment);

  private:
    /// Add visualization of a web segment, body based on primitives corresponding to the contact shapes.
    void AddWebVisualization(std::shared_ptr<ChBody> segment);

    std::vector<std::shared_ptr<ChBody>> m_web_segments;  ///< handles to track shoe's web segment bodies
    double m_seg_length;                                  ///< length of a web segment
    double m_seg_mass;                                    ///< mass of a web segment
    ChVector<> m_seg_inertia;                             ///< moments of inertia of a web segment

    friend class ChTrackAssemblyBandBushing;
};

/// Vector of handles to continuous band bushing-based track shoe subsystems.
typedef std::vector<std::shared_ptr<ChTrackShoeBandBushing>> ChTrackShoeBandBushingList;

/// @} vehicle_tracked_shoe

}  // end namespace vehicle
}  // end namespace chrono

#endif
