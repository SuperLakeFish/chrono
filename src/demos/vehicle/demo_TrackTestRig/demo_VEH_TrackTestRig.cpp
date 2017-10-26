// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
//
// =============================================================================

#include "chrono/core/ChFileutils.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_mkl/ChSolverMKL.h"

#include "chrono_vehicle/ChVehicleModelData.h"
#include "chrono_vehicle/utils/ChVehicleIrrApp.h"
#include "chrono_vehicle/tracked_vehicle/utils/ChTrackTestRig.h"
#include "chrono_vehicle/tracked_vehicle/utils/ChIrrGuiDriverTTR.h"

#include "chrono_vehicle/tracked_vehicle/track_assembly/TrackAssemblySinglePin.h"
////#include "chrono_vehicle/tracked_vehicle/track_assembly/TrackAssemblyDoublePin.h"
////#include "chrono_vehicle/tracked_vehicle/track_assembly/TrackAssemblyRigidCB.h"

#include "chrono_models/vehicle/m113/M113_TrackAssemblySinglePin.h"
#include "chrono_models/vehicle/m113/M113_TrackAssemblyDoublePin.h"
#include "chrono_models/vehicle/m113/M113_TrackAssemblyRigidCB.h"
#include "chrono_models/vehicle/m113/M113_TrackAssemblyRigidANCFCB.h"


#include "chrono_fea/ChElementShellANCF.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"


using namespace chrono;
using namespace chrono::vehicle;
using namespace chrono::vehicle::m113;
using namespace chrono::fea;

using std::cout;
using std::endl;

// =============================================================================
// USER SETTINGS
// =============================================================================

bool use_JSON = false;
std::string filename("M113/track_assembly/M113_TrackAssemblySinglePin_Left.json");

double post_limit = 0.2;

// Simulation step size
double step_size = 1e-4;

// Time interval between two render frames
double render_step_size = 1.0 / 500;

// Output (screenshot captures)
bool img_output = false;

const std::string out_dir = GetChronoOutputPath() + "TRACK_TEST_RIG";

// =============================================================================
int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    ChTrackTestRig* rig = nullptr;
    ChVector<> attach_loc(0, 1, 0);

    if (use_JSON) {
        rig = new ChTrackTestRig(vehicle::GetDataFile(filename), attach_loc);
    } else {
        VehicleSide side = LEFT;
        TrackShoeType type = TrackShoeType::RIGID_ANCF_CB;

        std::shared_ptr<ChTrackAssembly> track_assembly;
        switch (type) {
            case TrackShoeType::SINGLE_PIN: {
                auto assembly = std::make_shared<M113_TrackAssemblySinglePin>(side);
                track_assembly = assembly;
                break;
            }
            case TrackShoeType::DOUBLE_PIN: {
                auto assembly = std::make_shared<M113_TrackAssemblyDoublePin>(side);
                track_assembly = assembly;
                break;
            }
            case TrackShoeType::RIGID_CB: {
                auto assembly = std::make_shared<M113_TrackAssemblyRigidCB>(side);
                track_assembly = assembly;
                break;
            }
            case TrackShoeType::RIGID_ANCF_CB: {
                auto assembly = std::make_shared<M113_TrackAssemblyRigidANCFCB>(side);
                track_assembly = assembly;
                break;
            }
        }

        rig = new ChTrackTestRig(track_assembly, attach_loc, ChMaterialSurface::SMC);
    }

//-------------------------------------------------------------------
#if FALSE
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
    ChVector<> E(0.01e9, 0.01e9, 0.01e9);
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
    rig->GetSystem()->Add(my_mesh);

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
#endif

    //rig->GetSystem()->Set_G_acc(ChVector<>(0, 0, 0));
    rig->GetSystem()->SetSolverType(ChSolver::Type::SOR);
    rig->GetSystem()->SetMaxItersSolverSpeed(50);
    rig->GetSystem()->SetMaxItersSolverStab(50);
    rig->GetSystem()->SetTol(0);
    rig->GetSystem()->SetMaxPenetrationRecoverySpeed(1.5);
    rig->GetSystem()->SetMinBounceSpeed(2.0);
    rig->GetSystem()->SetSolverOverrelaxationParam(0.8);
    rig->GetSystem()->SetSolverSharpnessParam(1.0);
	rig->SetMaxTorque(6000);

	// Mark completion of system construction
	rig->GetSystem()->SetupInitial();

	auto mkl_solver = std::make_shared<ChSolverMKL<>>();
	rig->GetSystem()->SetSolver(mkl_solver);
	mkl_solver->SetSparsityPatternLock(false);
	rig->GetSystem()->Update();

	// HHT
	rig->GetSystem()->SetTimestepperType(ChTimestepper::Type::HHT);
	auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(rig->GetSystem()->GetTimestepper());
	mystepper->SetAlpha(-0.2);
	mystepper->SetMaxiters(200);
	mystepper->SetAbsTolerances(1e-03);
	mystepper->SetMode(ChTimestepperHHT::ACCELERATION);
	mystepper->SetScaling(true);
	mystepper->SetVerbose(true);
	mystepper->SetStepControl(true);

    

    ChVector<> rig_loc(0, 0, 2);
    ChQuaternion<> rig_rot(1, 0, 0, 0);
    rig->Initialize(ChCoordsys<>(rig_loc, rig_rot));

    rig->GetTrackAssembly()->SetSprocketVisualizationType(VisualizationType::PRIMITIVES);
    rig->GetTrackAssembly()->SetIdlerVisualizationType(VisualizationType::PRIMITIVES);
    rig->GetTrackAssembly()->SetRoadWheelAssemblyVisualizationType(VisualizationType::PRIMITIVES);
    rig->GetTrackAssembly()->SetRoadWheelVisualizationType(VisualizationType::PRIMITIVES);
    rig->GetTrackAssembly()->SetTrackShoeVisualizationType(VisualizationType::PRIMITIVES);

    ////rig->SetCollide(TrackedCollisionFlag::NONE);
    ////rig->SetCollide(TrackedCollisionFlag::SPROCKET_LEFT | TrackedCollisionFlag::SHOES_LEFT);
    ////rig->GetTrackAssembly()->GetSprocket()->GetGearBody()->SetCollide(false);

    // Create the vehicle Irrlicht application.
    ////ChVector<> target_point = rig->GetPostPosition();
    ////ChVector<> target_point = rig->GetTrackAssembly()->GetIdler()->GetWheelBody()->GetPos();
    ChVector<> target_point = rig->GetTrackAssembly()->GetSprocket()->GetGearBody()->GetPos();

    ChVehicleIrrApp app(rig, NULL, L"Suspension Test Rig");
    app.SetSkyBox();
    app.AddTypicalLights(irr::core::vector3df(30.f, -30.f, 100.f), irr::core::vector3df(30.f, 50.f, 100.f), 250, 130);
    app.SetChaseCamera(ChVector<>(0), 3.0, 0.0);
    app.SetChaseCameraPosition(target_point + ChVector<>(0, 3, 0));
    app.SetChaseCameraMultipliers(1e-4, 10);
    app.SetTimestep(step_size);
    app.AssetBindAll();
    app.AssetUpdateAll();

    // Create the driver system and set the time response for keyboard inputs.
    ChIrrGuiDriverTTR driver(app, post_limit);
    double steering_time = 1.0;      // time to go from 0 to max
    double displacement_time = 2.0;  // time to go from 0 to max applied post motion
    driver.SetSteeringDelta(render_step_size / steering_time);
    driver.SetDisplacementDelta(render_step_size / displacement_time * post_limit);
    driver.Initialize();

    // Initialize output
    if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
        std::cout << "Error creating directory " << out_dir << std::endl;
        return 1;
    }

    // ---------------
    // Simulation loop
    // ---------------

    // Inter-module communication data
    TerrainForces shoe_forces(1);

    // Number of simulation steps between two 3D view render frames
    int render_steps = (int)std::ceil(render_step_size / step_size);

    // Initialize simulation frame counter
    int step_number = 0;
    int render_frame = 0;

    while (app.GetDevice()->run()) {
        // Debugging output
        const ChFrameMoving<>& c_ref = rig->GetChassisBody()->GetFrame_REF_to_abs();
        const ChVector<>& i_pos_abs = rig->GetTrackAssembly()->GetIdler()->GetWheelBody()->GetPos();
        const ChVector<>& s_pos_abs = rig->GetTrackAssembly()->GetSprocket()->GetGearBody()->GetPos();
        ChVector<> i_pos_rel = c_ref.TransformPointParentToLocal(i_pos_abs);
        ChVector<> s_pos_rel = c_ref.TransformPointParentToLocal(s_pos_abs);
        ////cout << "Time: " << rig->GetSystem()->GetChTime() << endl;
        ////cout << "      idler:    " << i_pos_rel.x << "  " << i_pos_rel.y << "  " << i_pos_rel.z << endl;
        ////cout << "      sprocket: " << s_pos_rel.x << "  " << s_pos_rel.y << "  " << s_pos_rel.z << endl;

        // Render scene
        if (step_number % render_steps == 0) {
            app.BeginScene(true, true, irr::video::SColor(255, 140, 161, 192));
            app.DrawAll();
            app.EndScene();

            if (img_output && step_number > 1000) {
                char filename[100];
                sprintf(filename, "%s/img_%03d.jpg", out_dir.c_str(), render_frame + 1);
                app.WriteImageToFile(filename);
            }

            render_frame++;
        }

        // Collect output data from modules
        double throttle_input = driver.GetThrottle();
        double post_input = driver.GetDisplacement();

        // Update modules (process inputs from other modules)
        double time = rig->GetChTime();
        driver.Synchronize(time);
        rig->Synchronize(time, post_input, throttle_input, shoe_forces);
        app.Synchronize("", 0, throttle_input, 0);

        // Advance simulation for one timestep for all modules
        driver.Advance(step_size);
        rig->Advance(step_size);
        app.Advance(step_size);

        // Increment frame number
        step_number++;

		std::cout << "Step: " << step_number << "   Time: " << time << std::endl;
    }

    delete rig;

    return 0;
}
