//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

//
//   Demo code (advanced), about
//
//     - loading an Abaqus tetahedrom mesh
//     - apply a load to the mesh using an external tool, 
//       say CFD or SPH (here simulated as a function in this .cpp file)
//       that is perform a cosimulation.


#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChSystemDEM.h"
#include "chrono/physics/ChLoaderUV.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/lcp/ChLcpIterativeMINRES.h"

#include "chrono_fea/ChElementTetra_4.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChLoadContactSurfaceMesh.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono_irrlicht/ChIrrApp.h"


using namespace chrono;
using namespace fea;
using namespace irr;


// This function simulates the effect of an external program that 
// gets a triangle mesh and outputs the forces acting on the nodes.
// In a real cosimulation scenario, this procedure could even reside
// on a different computing node and manage inputs/outputs via MPI or such.

void PerformEsternalCosimulation(   const std::vector<ChVector<>>& input_vert_pos,
                                    const std::vector<ChVector<>>& input_vert_vel,
                                    const std::vector<ChVector<int>>& input_triangles,
                                    std::vector<ChVector<>>& vert_output_forces,
                                    std::vector<int>& vert_output_indexes
                                    ) {
    double ky = 10000; // upward stiffness
    double ry = 10; // upward damping
    // simple example: scan through all vertexes in the mesh, see if they sink below zero,
    // apply a penalty upward spring force if so.
    for (int iv = 0; iv < input_vert_pos.size(); ++iv) {
        if (input_vert_pos[iv].y < 0) {
            double yforce = - ky * input_vert_pos[iv].y - ry * input_vert_vel[iv].y;
            if (yforce > 0) {
                vert_output_forces.push_back(ChVector<>(0,yforce,0));
                vert_output_indexes.push_back(iv);
            }
        }
    }
    // note that: 
    // - we avoided addingfroces to  vert_output_forces  when force was zero.
    // - vert_output_forces has the same size of vert_output_indexes, maybe smaller than input_vert_pos 
}


int main(int argc, char* argv[]) {

    // Global parameter for tire:
    double tire_rad = 0.8;
    double tire_vel_z0 = -3;
    ChVector<> tire_center(0, 0.02+tire_rad, 0.5);
    ChMatrix33<> tire_alignment(Q_from_AngAxis(CH_C_PI, VECT_Y)); // create rotated 180� on y

    double tire_w0 = tire_vel_z0/tire_rad;


    // Create a Chrono::Engine physical system
    ChSystemDEM my_system;

    // Create the Irrlicht visualization (open the Irrlicht device,
    // bind a simple user interface, etc. etc.)
    ChIrrApp application(&my_system, L"FEA contacts", core::dimension2d<u32>(1280, 720), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(1, (f32)1.4, -1.2), core::vector3df(0, tire_rad, 0));
    application.AddLightWithShadow(core::vector3df(1.5, 5.5, -2.5), core::vector3df(0, 0, 0), 3, 2.2, 7.2, 40, 512,
                                   video::SColorf(0.8, 0.8, 1));

    //
    // CREATE THE PHYSICAL SYSTEM
    //

    // Create the surface material, containing information
    // about friction etc.

    ChSharedPtr<ChMaterialSurfaceDEM> mysurfmaterial (new ChMaterialSurfaceDEM);
    mysurfmaterial->SetKn(2e6);
    mysurfmaterial->SetKt(2e6);
    mysurfmaterial->SetGn(4200);
    mysurfmaterial->SetGt(4200);

  
    // Create a mesh, that is a container for groups
    // of FEA elements and their referenced nodes.

    ChSharedPtr<ChMesh> my_mesh(new ChMesh);
    my_system.Add(my_mesh); 


    // Create a material, that must be assigned to each solid element in the mesh,
    // and set its parameters

    ChSharedPtr<ChContinuumElastic> mmaterial(new ChContinuumElastic);
    mmaterial->Set_E(0.003e9);  // rubber 0.01e9, steel 200e9
    mmaterial->Set_v(0.4);
    mmaterial->Set_RayleighDampingK(0.004);
    mmaterial->Set_density(1000);


    // Load an ABAQUS .INP tetahedron mesh file from disk, defining a tetahedron mesh.
    // Note that not all features of INP files are supported. Also, quadratic tetahedrons are promoted to linear.
    // This is much easier than creating all nodes and elements via C++ programming.
    // Ex. you can generate these .INP files using Abaqus or exporting from the SolidWorks simulation tool.

    std::vector<std::vector<ChSharedPtr<ChNodeFEAbase> > > node_sets;

    try {
        my_mesh->LoadFromAbaqusFile(GetChronoDataFile("fea/tractor_wheel_coarse.INP").c_str(), mmaterial, node_sets, tire_center, tire_alignment);
    } catch (ChException myerr) {
        GetLog() << myerr.what();
        return 0;
    }


    // Create the contact surface(s). 
    // Use the AddFacesFromBoundary() to select automatically the outer skin of the tetrahedron mesh:

    ChSharedPtr<ChContactSurfaceMesh> mcontactsurf (new ChContactSurfaceMesh);
    my_mesh->AddContactSurface(mcontactsurf);
    
    mcontactsurf->AddFacesFromBoundary();

    mcontactsurf->SetMaterialSurface(mysurfmaterial); // by the way it is not needed because contacts will be emulated by cosimulation



    /// Create a mesh load for cosimulation, acting on the contact surface above
    /// (forces on nodes will be computed by an external procedure)

    ChSharedPtr<ChLoadContainer> mloadcontainer(new ChLoadContainer);
    my_system.Add(mloadcontainer);

    ChSharedPtr<ChLoadContactSurfaceMesh> mmeshload (new ChLoadContactSurfaceMesh(mcontactsurf));
    mloadcontainer->Add(mmeshload);
    


    //
    // Optional...  visualization
    //

    // ==Asset== attach a visualization of the FEM mesh.
    // This will automatically update a triangle mesh (a ChTriangleMeshShape
    // asset that is internally managed) by setting  proper
    // coordinates and vertex colours as in the FEM elements.

    ChSharedPtr<ChVisualizationFEAmesh> mvisualizemesh(new ChVisualizationFEAmesh(*(my_mesh.get_ptr())));
    mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
    mvisualizemesh->SetColorscaleMinMax(0.0, 10);
    mvisualizemesh->SetSmoothFaces(true);
    my_mesh->AddAsset(mvisualizemesh); 
    
 

    // ==IMPORTANT!== Use this function for adding a ChIrrNodeAsset to all items
    application.AssetBindAll();

    // ==IMPORTANT!== Use this function for 'converting' into Irrlicht meshes the assets
    application.AssetUpdateAll();

    // Use shadows in realtime view
    application.AddShadowAll();


    // ==IMPORTANT!== Mark completion of system construction
    my_system.SetupInitial();



    //
    // THE SOFT-REAL-TIME CYCLE
    //


    
        // Change solver to embedded MINRES
        // NOTE! it is strongly advised that you compile the optional MKL module 
        // if you need higher precision, and switch to its MKL solver - see demos for FEA & MKL.
    my_system.SetLcpSolverType(ChSystem::LCP_ITERATIVE_MINRES);     
    my_system.SetIterLCPwarmStarting(true);  // this helps a lot to speedup convergence in this class of problems
    my_system.SetIterLCPmaxItersSpeed(40);
    my_system.SetTolForce(1e-10);  
    

    // Change type of integrator:
    my_system.SetIntegrationType(chrono::ChSystem::INT_EULER_IMPLICIT_LINEARIZED);  // fast, less precise


    application.SetTimestep(0.005);

    while (application.GetDevice()->run()) {
        application.BeginScene();

        application.DrawAll();

        application.DoStep();

        // -------------------------------------------------------------------------
        // Here do the cosimulation
        // A <--> B  
        // For example, A is this main program, and B can be an external
        // program, ex. a CFD or SPH simulation tool. 
        // The idea is that A --> B communicates the mesh position,
        // then A <-- B receives the computed forces to be applied at nodes.
        // In this example, to keep things simple, B is just a simple C function
        // in this .cpp file.

        std::vector<ChVector<>> vert_pos;
        std::vector<ChVector<>> vert_vel;
        std::vector<ChVector<int>> triangles;
        std::vector<ChVector<>> vert_forces;
        std::vector<int> vert_indexes;

        mmeshload->OutputSimpleMesh(vert_pos, vert_vel, triangles);

        PerformEsternalCosimulation(vert_pos, 
                                    vert_vel, 
                                    triangles, 
                                    vert_forces, 
                                    vert_indexes);

        mmeshload->InputSimpleForces(vert_forces, vert_indexes);

        // End of cosimulation block
        // -------------------------------------------------------------------------


        // now, just for debugging and some fun, draw some triangles
        // (only those that have a vertex that has a force applied):
        for (int it= 0;it < triangles.size(); ++it) {
            bool vert_hit = false;
            for (int io = 0; io < vert_indexes.size(); ++io) {
                if (triangles[it].x == vert_indexes[io] || triangles[it].y == vert_indexes[io] || triangles[it].z == vert_indexes[io])
                    vert_hit = true;
            }
            if (vert_hit == true) {
                std::vector<chrono::ChVector<> > fourpoints = { vert_pos[triangles[it].x], 
                                                                vert_pos[triangles[it].y], 
                                                                vert_pos[triangles[it].z],
                                                                vert_pos[triangles[it].x]};
                ChIrrTools::drawPolyline(application.GetVideoDriver(), fourpoints, irr::video::SColor(255,240,200,0), true); 
            }
        }


        ChIrrTools::drawGrid(application.GetVideoDriver(), 0.1, 0.1, 20, 20,
                                 ChCoordsys<>(VNULL, CH_C_PI_2, VECT_X), video::SColor(50, 90, 90, 90), true);

        application.EndScene();
    }

    return 0;
}