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
// Authors: Michael Taylor
// =============================================================================
//
// Demo on using ANCF shell elements
//
// =============================================================================


#include "chrono_mkl/ChSolverMKL.h"
#include "chrono/solver/ChSolverMINRES.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono_fea/ChElementShellANCF.h"
#include "chrono_fea/ChElementShellANCF_8.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono_irrlicht/ChIrrApp.h"

using namespace chrono;
using namespace chrono::fea;
using namespace chrono::irrlicht;

using namespace irr;

int main(int argc, char* argv[]) {
	GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

	double time_step = 1e-3;

	ChSystemSMC my_system;
	my_system.Set_G_acc(ChVector<>(0, 0, -9.8*0.0));

	// Create the Irrlicht visualization (open the Irrlicht device, bind a simple user interface, etc.)
	ChIrrApp application(&my_system, L"ANCF Shells", core::dimension2d<u32>(800, 600), false, true);

	// Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
	application.AddTypicalLogo();
	application.AddTypicalSky();
	application.AddTypicalLights();
	application.AddTypicalCamera(core::vector3df(-0.4f, -0.3f, 0.0f),  // camera location
		core::vector3df(0.0f, 0.5f, -0.1f));  // "look at" location

	GetLog() << "-----------------------------------------------------------\n";
	GetLog() << "-----------------------------------------------------------\n";
	GetLog() << "  ANCF 8 node Shell Elements demo for CB pre-stressed use \n";
	GetLog() << "-----------------------------------------------------------\n";

	//Problem Geometry
	double web_angle = -30 * CH_C_DEG_TO_RAD;
	double length = 24 * 25.4/1000;
	double width = 12 * 25.4 / 1000;
	double thickness = 1 * 25.4 / 1000;
	
	int num_elements_length = 10;
	int num_elements_width = 10;
	int num_elements_thickness = 1;

    int N_x_edge = 2 * num_elements_length + 1;
    int N_y_edge = 2 * num_elements_width + 1;
    int N_x_mid = num_elements_length + 1;
    int N_y_mid = num_elements_width + 1;

	double dx = length * std::cos(web_angle) / (2*num_elements_length);
	double dy = width / (2*num_elements_width);
	double dz = length * std::sin(web_angle) / (2*num_elements_length);

	// Create a mesh, that is a container for groups of elements and their referenced nodes.
	auto my_mesh = std::make_shared<ChMesh>();

	// Create and add the nodes
	for (int x_idx = 0; x_idx < N_x_edge; x_idx++) {
		for (int y_idx = 0; y_idx < N_y_edge; y_idx++) {
            if ((x_idx % 2 == 1) && (y_idx % 2 == 1))
                continue;

			// Node location
			double loc_x = x_idx * dx;
			double loc_y = y_idx * dy;
			double loc_z = x_idx * dz;

			// Node direction
			double dir_x = std::cos(web_angle + CH_C_PI_2);
			double dir_y = 0;
			double dir_z = std::sin(web_angle + CH_C_PI_2);

            // Node direction derivative
            double curvz_x = 0.0;
            double curvz_y = 0.0;
            double curvz_z = 0.0;

			// Create the node
			auto node = std::make_shared<ChNodeFEAxyzDD>(ChVector<>(loc_x, loc_y, loc_z), ChVector<>(dir_x, dir_y, dir_z),
                ChVector<>(curvz_x, curvz_y, curvz_z));

			node->SetMass(0);

			// Add node to mesh
			my_mesh->AddNode(node);
		}
	}


	// Create an orthotropic material.
	// All layers for all elements share the same material.
	double rho = 500;
	ChVector<> E(2.1e7, 2.1e7, 2.1e7);
	ChVector<> nu(0.3, 0.3, 0.3);
	ChVector<> G(8.0769231e6, 8.0769231e6, 8.0769231e6);
	auto mat = std::make_shared<ChMaterialShellANCF_8>(rho, E, nu, G);

	// Create the elements
	for (int x_idx = 0; x_idx < num_elements_length; x_idx++) {
		for (int y_idx = 0; y_idx < num_elements_width; y_idx++) {
			// Adjacent nodes
            /// The node numbering is in ccw fashion as in the following scheme:
            ///         v
            ///         ^
            /// D o-----G-----o C
            ///   |     |     |
            /// --H-----+-----F-> u
            ///   |     |     |
            /// A o-----E-----o B

			int node0 = 2*y_idx + x_idx * (N_y_edge + N_y_mid);
			int node1 = 2*y_idx + (x_idx + 1 ) * (N_y_edge + N_y_mid);
			int node2 = 2*(y_idx + 1) + (x_idx + 1) * (N_y_edge + N_y_mid);
			int node3 = 2*(y_idx + 1) + x_idx * (N_y_edge + N_y_mid);

            int node4 = 1 * y_idx + x_idx * (N_y_edge + N_y_mid) + N_y_edge;
            int node5 = 2 * y_idx + (x_idx + 1) * (N_y_edge + N_y_mid) + 1;
            int node6 = 1 * y_idx + x_idx * (N_y_edge + N_y_mid) + N_y_edge + 1;
            int node7 = 2 * y_idx + x_idx * (N_y_edge + N_y_mid) + 1;


			// Create the element and set its nodes.
            auto element = std::make_shared<ChElementShellANCF_8>();
            element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzDD>(my_mesh->GetNode(node0)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(my_mesh->GetNode(node1)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(my_mesh->GetNode(node2)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(my_mesh->GetNode(node3)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(my_mesh->GetNode(node4)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(my_mesh->GetNode(node5)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(my_mesh->GetNode(node6)),
                std::dynamic_pointer_cast<ChNodeFEAxyzDD>(my_mesh->GetNode(node7)));

			// Set element dimensions
			element->SetDimensions(2*std::sqrt(dx*dx + dz*dz), 2*dy);

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
	my_system.Add(my_mesh);

	auto fixed_body = std::make_shared<ChBodyEasyCylinder>(0.01, width*1.05, 500);
    fixed_body->SetPos(ChVector<>(0.0, (width*1.05)/2.0, 0.0));
    fixed_body->SetBodyFixed(true);
	my_system.Add(fixed_body);

	auto tip_body = std::make_shared<ChBodyEasyCylinder>(0.01, width*1.05, 500);
	tip_body->SetPos(ChVector<>(2*dx*num_elements_length, (width*1.05) / 2.0, 2*dz*num_elements_length));
    tip_body->SetBodyFixed(true);
    my_system.Add(tip_body);

    // Change the gradient on the boundary nodes that will connect to the first fixed body
    // and then connect those nodes to the body
    for (int y_idx = 0; y_idx < N_y_edge; y_idx++) {
        int node_idx = y_idx;
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzDD>(my_mesh->GetNode(node_idx));

        node->SetD(ChVector<>(0.0, 0.0, 1.0));

        auto constraintxyz = std::make_shared<ChLinkPointFrame>();
        constraintxyz->Initialize(node, fixed_body);
        my_system.Add(constraintxyz);

        auto constraintD = std::make_shared<ChLinkDirFrame>();
        constraintD->Initialize(node, fixed_body);
        my_system.Add(constraintD);
    }

    // Change the gradient on the boundary nodes that will connect to the second fixed body
    // and then connect those nodes to the body
    for (int y_idx = 0; y_idx < N_y_edge; y_idx++) {
        int node_idx = y_idx + num_elements_length * (N_y_edge + N_y_mid);
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzDD>(my_mesh->GetNode(node_idx));

        node->SetD(ChVector<>(std::sin(60 * CH_C_DEG_TO_RAD), 0.0, std::cos(60*CH_C_DEG_TO_RAD)));

        auto constraintxyz = std::make_shared<ChLinkPointFrame>();
        constraintxyz->Initialize(node, tip_body);
        my_system.Add(constraintxyz);

        auto constraintD = std::make_shared<ChLinkDirFrame>();
        constraintD->Initialize(node, tip_body);
        my_system.Add(constraintD);
    }


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

	application.AssetBindAll();
	application.AssetUpdateAll();

	// ----------------------------------
	// Perform a dynamic time integration
	// ----------------------------------

	// Mark completion of system construction
	my_system.SetupInitial();


	// Set up solver
	auto mkl_solver = std::make_shared<ChSolverMKL<>>();
	my_system.SetSolver(mkl_solver);
	mkl_solver->SetSparsityPatternLock(false);
	my_system.Update();

	// HHT
	my_system.SetTimestepperType(ChTimestepper::Type::HHT);
	auto mystepper = std::dynamic_pointer_cast<ChTimestepperHHT>(my_system.GetTimestepper());
	mystepper->SetAlpha(-0.2);
	mystepper->SetMaxiters(200);
	mystepper->SetAbsTolerances(1e-5);
	mystepper->SetMode(ChTimestepperHHT::ACCELERATION);
	mystepper->SetScaling(true);
	mystepper->SetVerbose(false);
	mystepper->SetStepControl(true);
	application.SetTimestep(time_step);

	size_t step_num = 0;
	while (application.GetDevice()->run()) {
		application.BeginScene();
		application.DrawAll();
		application.DoStep();
		application.EndScene();

		if ((step_num % 10)==0) {
			std::cout << "Time: " << my_system.GetChTime() << std::endl;
		}
		
		step_num++;
	}

	return 0;
}
