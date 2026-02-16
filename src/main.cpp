#include "mfem.hpp"
#include "dg_block_inverse.hpp"
#include "maxwell_evolution.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace mfem;
using namespace std;

int main(int argc, char *argv[]){
  const char *mesh_file = "task.mesh";
   int ref_levels = 0;
   int order = 1;
   double dt = 1e-13;
   double t_final = 2.5e-8;
   int vis_steps =10;
   bool visualization = true;
   bool paraview = true;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine", "Number of times to refine the mesh.");
   args.AddOption(&order, "-o", "--order", "Finite element order (polynomial degree).");
   args.AddOption(&dt, "-dt", "--time-step", "Time step.");
   args.AddOption(&t_final, "-tf", "--final-time", "Final time of simulation.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis", "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&paraview, "-paraview", "--paraview-datafiles", "-no-paraview",
                  "--no-paraview-datafiles",
                  "Save data files for ParaView (paraview.org) visualization.");
   args.AddOption(&vis_steps, "-vs", "--vis_steps", "Print visualizations every #N steps");
   args.Parse();
   if (!args.Good()) { args.PrintUsage(cout); return 1; }
   args.PrintOptions(cout);

   Mesh mesh(mesh_file, 1, 1);
   for (int i = 0; i < ref_levels; i++) { mesh.UniformRefinement(); }

   const int dim = mesh.Dimension();
   DG_FECollection fec(order, dim);
   FiniteElementSpace Hx(&mesh, &fec, 1, Ordering::byVDIM);
   FiniteElementSpace Hy(&mesh, &fec, 1, Ordering::byVDIM);
   FiniteElementSpace Ez(&mesh, &fec, 1, Ordering::byVDIM);
   cout << "Dofs = " << 3*Hx.GetVSize() << endl;
      cout << "Boundary attributes max = "
     << mesh.bdr_attributes.Max() << endl;

   GridFunction Hx_gf(&Hx);
   GridFunction Hy_gf(&Hy);
   GridFunction Ez_gf(&Ez);

   //Assign Mesh attributes for inlet outlet and pec
   int max_bdr = mesh.bdr_attributes.Max();

   Array<int> attr_inlet(max_bdr); attr_inlet = 0;
   if(max_bdr >= 1) attr_inlet[0] = 1;

   Array<int> attr_outlet(max_bdr); attr_outlet = 0;
   if(max_bdr >= 2) attr_outlet[1] = 1;

   Array<int> attr_walls(max_bdr); attr_walls = 0;
   if(max_bdr >= 3) attr_walls[2] = 1;

   double epsilon0 = 8.8541878176e-12;
   double mu0 = 1.2566370614e-6;
   double Z0 = sqrt(mu0 / epsilon0);

   double f0 = 11e9;
   double H_height = 0.02;
   double omega = 2.0 * M_PI * f0;
   double k0 = omega * sqrt(mu0 * epsilon0);
   double kc = M_PI / H_height;

   if (k0 <= kc) {
        cout << "Error: Frequency below cutoff!" << endl;
        return 1;
    }


   ConstantCoefficient one(1.0);
   ConstantCoefficient minusone(-1.0);
   ConstantCoefficient mu(mu0);
   ConstantCoefficient epsilon(epsilon0);
   ConstantCoefficient Z(1.0*Z0);
   ConstantCoefficient Y(1.0/(Z0));

   Vector par_x(2); par_x = 0.0; par_x(0) = 1.0; VectorConstantCoefficient dell_x(par_x);
   Vector par_y(2); par_y = 0.0; par_y(1) = 1.0; VectorConstantCoefficient dell_y(par_y);
   ConstantCoefficient alpha(0.5);
   ConstantCoefficient minusalpha(-0.5);

   //Bilinear Forms for the Mass Matrix
   BilinearForm mhx(&Hx);
   BilinearForm mhy(&Hy);
   BilinearForm mez(&Ez);

   DGBlockInverse inv_mhx(Hx, mu0);
   DGBlockInverse inv_mhy(Hy, mu0);
   DGBlockInverse inv_mez(Ez, epsilon0);

   cout << "Mass Matrix Assembled" << endl;

   //Bilinear Forms for the K and F matrices
   BilinearForm k_hx_ez(&Hx);
   BilinearForm k_hy_ez(&Hy);
   BilinearForm k_ez_hx(&Ez);
   BilinearForm k_ez_hy(&Ez);

   k_hx_ez.AddDomainIntegrator(new TransposeIntegrator(new DerivativeIntegrator(one,1)));
   k_hy_ez.AddDomainIntegrator(new TransposeIntegrator(new DerivativeIntegrator(minusone,0)));
   k_ez_hx.AddDomainIntegrator(new TransposeIntegrator(new DerivativeIntegrator(one,1)));
   k_ez_hy.AddDomainIntegrator(new TransposeIntegrator(new DerivativeIntegrator(minusone,0)));

   cout << "Directional Derivative added" << endl;

   k_hx_ez.AddInteriorFaceIntegrator(new DGTraceIntegrator(dell_y,-1.0,0.0));
   k_hy_ez.AddInteriorFaceIntegrator(new DGTraceIntegrator(dell_x,1.0,0.0));
   k_ez_hx.AddInteriorFaceIntegrator(new DGTraceIntegrator(dell_y,-1.0, 0.0));
   k_ez_hy.AddInteriorFaceIntegrator(new DGTraceIntegrator(dell_x,1.0, 0.0));

   cout << "DG Central Flux added"<< endl;

   k_ez_hx.AddBdrFaceIntegrator(new DGTraceIntegrator(dell_y,-2.0, 0.0), attr_walls);
   k_ez_hy.AddBdrFaceIntegrator(new DGTraceIntegrator(dell_x,2.0, 0.0), attr_walls);

   cout << "PEC added" << endl;

   k_hx_ez.Assemble(); k_hx_ez.Finalize();
   k_hy_ez.Assemble(); k_hy_ez.Finalize();
   k_ez_hx.Assemble(); k_ez_hx.Finalize();
   k_ez_hy.Assemble(); k_ez_hy.Finalize();

   cout << "Assembled Flux and Directional bilinear form" << endl;

   //Bilinear Forms for the ABC
   BilinearForm abc_ez_ez(&Ez);
   BilinearForm abc_hy_hy(&Hy);

   abc_ez_ez.AddBdrFaceIntegrator(new BoundaryMassIntegrator(Y), attr_outlet);
   abc_ez_ez.AddBdrFaceIntegrator(new BoundaryMassIntegrator(Y), attr_inlet);

   abc_hy_hy.AddBdrFaceIntegrator(new BoundaryMassIntegrator(Z), attr_outlet);
   abc_hy_hy.AddBdrFaceIntegrator(new BoundaryMassIntegrator(Z), attr_inlet);

   abc_ez_ez.Assemble(); abc_ez_ez.Finalize();
   abc_hy_hy.Assemble(); abc_hy_hy.Finalize();

   cout << "Assembled ABC" << endl;

   //Create Big Block Matrix (F+K)
   Array<int> block_offsets(4);
   block_offsets[0] = 0;
   block_offsets[1] = Hx.GetVSize();
   block_offsets[2] = block_offsets[1] + Hy.GetVSize();
   block_offsets[3] = block_offsets[2] + Ez.GetVSize();

   BlockOperator A(block_offsets);

   A.SetBlock(0, 2, &k_hx_ez.SpMat());
   A.SetBlock(1, 2, &k_hy_ez.SpMat());
   A.SetBlock(2, 0, &k_ez_hx.SpMat());
   A.SetBlock(2, 1, &k_ez_hy.SpMat());
   A.SetBlock(2, 2, &abc_ez_ez.SpMat());
   A.SetBlock(1, 1, &abc_hy_hy.SpMat());

   //Solution Vector
   Vector Hx_vec(Hx.GetVSize()); Hx_vec = 0.0;
   Vector Hy_vec(Hy.GetVSize()); Hy_vec = 0.0;
   Vector Ez_vec(Ez.GetVSize()); Ez_vec = 0.0;

   cout << "Setup Big A Matrix" << endl;

   //Set boundary to zero, not necessary really
   ConstantCoefficient zero_coeff(0.0);
   Ez_gf.ProjectCoefficient(zero_coeff, attr_walls);
   Ez_gf.GetTrueDofs(Ez_vec);

   // Create BlockVector for full system
   BlockVector U(block_offsets);
   U.GetBlock(0) = Hx_vec;
   U.GetBlock(1) = Hy_vec;
   U.GetBlock(2) = Ez_vec;

   cout << "Setup Big U Vector" << endl;

   // Define Mode of the incoming wave
   FunctionCoefficient inlet_coeff([](const Vector &x){
     return sin(M_PI * (x(1)) / 0.02);
   });

   //E_inc linear form. Are later multiplied by the value in the TimeDependent Operator
   LinearForm f_inlet_ez(&Ez);
   f_inlet_ez.AddBdrFaceIntegrator(new BoundaryLFIntegrator(inlet_coeff), attr_inlet);
   f_inlet_ez.Assemble();
   LinearForm f_inlet_hy(&Hy);
   f_inlet_hy.AddBdrFaceIntegrator(new BoundaryLFIntegrator(inlet_coeff), attr_inlet);
   f_inlet_hy.Assemble();

   cout << "Assembled Linear Forms for the Inlet" << endl;

   // Time Dependent Operator
   MaxwellEvolution maxwell(A, inv_mhx, inv_mhy, inv_mez, &f_inlet_ez, &f_inlet_hy, block_offsets);

   double t = 0.0;
   int step = 0;
   int nsteps = int(t_final / dt);

   RK4Solver ode_solver;
   ode_solver.Init(maxwell);

   ostringstream sol_name;
   sol_name << "Ez";

   //Set Up Data collection for output
   ParaViewDataCollection paraview_dc("MaxwellSim", &mesh);
   paraview_dc.SetLevelsOfDetail(order);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.RegisterField("Hx", &Hx_gf);
   paraview_dc.RegisterField("Hy", &Hy_gf);
   paraview_dc.RegisterField("Ez", &Ez_gf);
   paraview_dc.SetCycle(0);
   paraview_dc.SetTime(0.0);

   // Save initial condition
   paraview_dc.SetCycle(0);
   paraview_dc.SetCycle(0.0);
   paraview_dc.Save();

   cout << "Sucessfully entered the loop" << endl;

   for (step = 1; step < nsteps; step++)
     {
       ode_solver.Step(U,t,dt);

        if (step % vis_steps == 0)
        {

            Hx_gf = U.GetBlock(0);
            Hy_gf = U.GetBlock(1);
            Ez_gf = U.GetBlock(2);

            paraview_dc.SetCycle(step);
            paraview_dc.SetTime(t);
            paraview_dc.Save();

            cout << "Step " << step << ", t = " << t << endl;

        }
    }
    cout << "Exited Loop" << endl;
    return 0;
}
