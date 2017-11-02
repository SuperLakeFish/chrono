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
	GetLog() << "     ANCF Shell Elements demo for CB pre-stressed use \n";
	GetLog() << "-----------------------------------------------------------\n";

	//Problem Geometry
    double radius = 8 * 25.4/1000;
    double thickness = 1.75 * 25.4 / 1000;
    double width = 12.5 * 25.4/1000;
    double box_length = 3.0 * 25.4 / 1000;
    int num_sections = 8;
	
	int num_elements_length = 3;
	int num_elements_width = 4;

	int N_x = num_elements_length + 1;
	int N_y = num_elements_width + 1;

    //Other inferred geometry
    double section_angle = CH_C_2PI / num_sections;
    double box_half_angle = std::asin((box_length / 2.0) / radius);
    double box_center_rad = radius*std::sin(box_half_angle)/std::tan(box_half_angle);
    double shell_angle = section_angle - 2.0 * box_half_angle;
    double shell_length = 2.0*radius*std::sin(shell_angle / 2.0);

	double dx = shell_length / num_elements_length;
	double dy = width / num_elements_width;
	double dz = thickness;

	// Number of elements in the z direction is considered as 1
	int TotalNumElementsPerSection = num_elements_length * num_elements_width;
	int TotalNumNodesPerSection = N_x * N_y;

    //Create an orthotropic material.
    //All layers for all elements share the same material.
    double rho = 500;
    ChVector<> E(2.1e7, 2.1e7, 2.1e7);
    ChVector<> nu(0.3, 0.3, 0.3);
    ChVector<> G(8.0769231e6, 8.0769231e6, 8.0769231e6);
    auto mat = std::make_shared<ChMaterialShellANCF>(rho, E, nu, G);

    std::vector<std::shared_ptr<ChBodyEasyBox>> boxes;
    std::vector<std::shared_ptr<ChMesh>> meshes;

    //Create the segments
    for (int seg_idx = 0; seg_idx < num_sections; seg_idx++) {
        
        //Create the rigid box
        boxes.push_back(std::make_shared<ChBodyEasyBox>(box_length, width, thickness, 500));
        boxes[seg_idx]->SetPos(box_center_rad * ChVector<>(std::cos(seg_idx*section_angle), 0.0, std::sin(seg_idx*section_angle)));
        boxes[seg_idx]->SetRot(Q_from_AngY(-CH_C_PI_2 - seg_idx*section_angle));
        if (seg_idx == 0) {
            boxes[seg_idx]->SetBodyFixed(true);
        }
        my_system.Add(boxes[seg_idx]);
        
        // Create a mesh, that is a container for groups of elements and their referenced nodes.
        meshes.push_back(std::make_shared<ChMesh>());
        
        ChVector<> xdir(std::cos((0.5 + seg_idx)*section_angle + CH_C_PI_2), 0.0, std::sin((0.5 + seg_idx)*section_angle + CH_C_PI_2));
        ChVector<> ydir(0.0, 1.0, 0.0);
        ChVector<> zdir(std::cos((0.5 + seg_idx)*section_angle + CH_C_PI), 0.0, std::sin((0.5 + seg_idx)*section_angle + CH_C_PI));

        ChVector<> StartingPoint(radius*std::cos(box_half_angle + seg_idx*section_angle), -0.5*width, radius*std::sin(box_half_angle + seg_idx*section_angle));

        // Create and add the nodes
        for (int x_idx = 0; x_idx < N_x; x_idx++) {
        	for (int y_idx = 0; y_idx < N_y; y_idx++) {
        	
        		// Node location
                ChVector<> node_loc = StartingPoint + (xdir * dx * x_idx) + (ydir * dy * y_idx);

        		// Node direction
                ChVector<> node_dir = zdir;

        		// Create the node
        		auto node = std::make_shared<ChNodeFEAxyzD>(node_loc, node_dir);

        		node->SetMass(0);

        		// Add node to mesh
                meshes[seg_idx]->AddNode(node);
        	}
        }


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
        		element->SetNodes(std::dynamic_pointer_cast<ChNodeFEAxyzD>(meshes[seg_idx]->GetNode(node0)),
        			std::dynamic_pointer_cast<ChNodeFEAxyzD>(meshes[seg_idx]->GetNode(node1)),
        			std::dynamic_pointer_cast<ChNodeFEAxyzD>(meshes[seg_idx]->GetNode(node2)),
        			std::dynamic_pointer_cast<ChNodeFEAxyzD>(meshes[seg_idx]->GetNode(node3)));

        		// Set element dimensions
        		element->SetDimensions(dx, dy);

        		// Add a single layers with a fiber angle of 0 degrees.
        		element->AddLayer(dz, 0 * CH_C_DEG_TO_RAD, mat);

        		// Set other element properties
        		element->SetAlphaDamp(0.05);    // Structural damping for this element
        		element->SetGravityOn(false);  // turn internal gravitational force calculation off

        		    // Add element to mesh
                meshes[seg_idx]->AddElement(element);
        	}
        }
        
        // -------------------------------------
        // Options for visualization in irrlicht
        // -------------------------------------

        auto mvisualizemesh = std::make_shared<ChVisualizationFEAmesh>(*(meshes[seg_idx].get()));
        mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
        mvisualizemesh->SetColorscaleMinMax(0.0, 5.50);
        mvisualizemesh->SetShrinkElements(true, 0.85);
        mvisualizemesh->SetSmoothFaces(true);
        meshes[seg_idx]->AddAsset(mvisualizemesh);

        auto mvisualizemeshref = std::make_shared<ChVisualizationFEAmesh>(*(meshes[seg_idx].get()));
        mvisualizemeshref->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
        mvisualizemeshref->SetWireframe(true);
        mvisualizemeshref->SetDrawInUndeformedReference(true);
        meshes[seg_idx]->AddAsset(mvisualizemeshref);

        auto mvisualizemeshC = std::make_shared<ChVisualizationFEAmesh>(*(meshes[seg_idx].get()));
        mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);
        mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
        mvisualizemeshC->SetSymbolsThickness(0.004);
        meshes[seg_idx]->AddAsset(mvisualizemeshC);

        auto mvisualizemeshD = std::make_shared<ChVisualizationFEAmesh>(*(meshes[seg_idx].get()));
        // mvisualizemeshD->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_VECT_SPEED);
        mvisualizemeshD->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_ELEM_TENS_STRAIN);
        mvisualizemeshD->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
        mvisualizemeshD->SetSymbolsScale(1);
        mvisualizemeshD->SetColorscaleMinMax(-0.5, 5);
        mvisualizemeshD->SetZbufferHide(false);
        meshes[seg_idx]->AddAsset(mvisualizemeshD);


        // Add the mesh to the system
        my_system.Add(meshes[seg_idx]);
    }

    //Link the segments
    for (int seg_idx = 0; seg_idx < num_sections; seg_idx++) {

        int next_idx = seg_idx + 1;
        if (next_idx == num_sections) {
            next_idx = 0;
        }

        auto rot = boxes[seg_idx]->GetRot();
        ChVector<> zdir = rot.GetZaxis();
        auto rot_next = boxes[next_idx]->GetRot();
        ChVector<> zdir_next = rot_next.GetZaxis();

        // Change the gradient on the boundary nodes that will connect to the first fixed body
        // and then connect those nodes to the body
        for (int y_idx = 0; y_idx < N_y; y_idx++) {
            int node_idx = y_idx + 0 * N_y;
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyzD>(meshes[seg_idx]->GetNode(node_idx));

            node->SetD(zdir);

            auto constraintxyz = std::make_shared<ChLinkPointFrame>();
            constraintxyz->Initialize(node, boxes[seg_idx]);
            my_system.Add(constraintxyz);

            auto constraintD = std::make_shared<ChLinkDirFrame>();
            constraintD->Initialize(node, boxes[seg_idx]);
            my_system.Add(constraintD);
        }

        // Change the gradient on the boundary nodes that will connect to the second fixed body
        // and then connect those nodes to the body
        for (int y_idx = 0; y_idx < N_y; y_idx++) {
            int node_idx = y_idx + num_elements_length * N_y;
            auto node = std::dynamic_pointer_cast<ChNodeFEAxyzD>(meshes[seg_idx]->GetNode(node_idx));

            node->SetD(zdir_next);

            auto constraintxyz = std::make_shared<ChLinkPointFrame>();
            constraintxyz->Initialize(node, boxes[next_idx]);
            my_system.Add(constraintxyz);

            auto constraintD = std::make_shared<ChLinkDirFrame>();
            constraintD->Initialize(node, boxes[next_idx]);
            my_system.Add(constraintD);
        }

    }

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
