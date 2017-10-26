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
// Authors: Milad Rakhsha, Radu Serban
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
	GetLog() << "     ANCF Shell Elements demo for CB pre-stressed use \n";
	GetLog() << "-----------------------------------------------------------\n";

	//Problem Geometry
	double web_angle = -30 * CH_C_DEG_TO_RAD;
	double length = 24 * 25.4/1000;
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


	//// Create an orthotropic material.
	//// All layers for all elements share the same material.
	//double rho = 500;
	//ChVector<> E(2.1e7, 2.1e7, 2.1e7);
	//ChVector<> nu(0.3, 0.3, 0.3);
	//ChVector<> G(8.0769231e6, 8.0769231e6, 8.0769231e6);
	//auto mat = std::make_shared<ChMaterialShellANCF>(rho, E, nu, G);

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
			int node1 = y_idx + (x_idx + 1 ) * N_y;
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
	my_system.Add(my_mesh);

	auto fixed_body = std::make_shared<ChBodyEasyCylinder>(0.01, width*1.05, 500);
    fixed_body->SetPos(ChVector<>(0.0, (width*1.05)/2.0, 0.0));
    fixed_body->SetBodyFixed(true);
	my_system.Add(fixed_body);

	auto tip_body = std::make_shared<ChBodyEasyCylinder>(0.01, width*1.05, 500);
	tip_body->SetPos(ChVector<>(dx*num_elements_length, (width*1.05) / 2.0, dz*num_elements_length));
    tip_body->SetBodyFixed(true);
    my_system.Add(tip_body);

    // Change the gradient on the boundary nodes that will connect to the first fixed body
    // and then connect those nodes to the body
    for (int y_idx = 0; y_idx < num_elements_width; y_idx++) {
        int node_idx = y_idx + 0 * N_y;
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node_idx));

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
    for (int y_idx = 0; y_idx < num_elements_width; y_idx++) {
        int node_idx = y_idx + num_elements_length * N_y;
        auto node = std::dynamic_pointer_cast<ChNodeFEAxyzD>(my_mesh->GetNode(node_idx));

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
