#include <gmsh.h>
#include <vector>
#include <cassert>
#include <string>
#include <iostream>

int main(int argc, char **argv){
  gmsh::initialize();
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  gmsh::option::setNumber("Mesh.Algorithm", 9);

  gmsh::model::add("waveguide");
  std::string name = "waveguide_mesh.msh";

  double lc = 0.005;
  double L_in = 0.1;
  double L_out = 0.079;
  double H = 0.02;
  double gap_width = 0.003;
  double gap_height = 0.005;
  double iris_width = 0.015;
  double scale = 3.0;

  double L_total = 2 * gap_width + L_out + L_in + iris_width;

  for (int i =1; i < argc; ++i) {
    std::string arg = argv[i];
    if ((arg == "-lc" || arg == "-mesh") && i + 1 < argc) {
      lc = std::atof(argv[++i]);
    } else if(arg == "-lin" &&  i + 1 < argc) {
      L_in = std::atof(argv[++i]);
    } else if(arg == "-lout" &&  i + 1 < argc) {
      L_out = std::atof(argv[++i]);
    } else if(arg == "-H" && i + 1 < argc) {
      H = std::atof(argv[++i]);
    } else if(arg == "-gap_w" &&  i + 1 < argc) {
      gap_width = std::atof(argv[++i]);
    } else if(arg == "-gap_h" &&  i + 1 < argc) {
      gap_height = std::atof(argv[++i]);
    } else if(arg == "-iris_w" &&  i + 1 < argc) {
      iris_width = std::atof(argv[++i]);
    } else if (arg == "-name" && i+1<argc) {
      name = argv[++i];
     } else if (arg == "-scale" && i+1<argc) {
        scale = std::atof(argv[++i]);
    }else if (arg == "-help") {
      std::cout << "Usage: " << argv[0] << " [options]\n"
               << "Options:\n"
               << "  -lc <val>       Set mesh size (default 0.005)\n"
               << "  -lin <val>      Set input waveguide size (default 0.1)\n"
               << "  -lout <val>     Set output waveguide size (default 0.079)\n"
               << "  -H <val>        Set waveguide height (default 0.02)\n"
               << "  -gap_w <val>    Set iris input seperator width (default 0.003)\n"
               << "  -gap_h <val>    Set iris input seperator height (default 0.005)\n"
               << "  -iris_w <val>   Set iris width (default 0.015)\n"
               <<"  -name <string>  Set Name to be saved to (Ending in .msh, default waveguide_mesh.msh)\n"
               <<" -scale <val>     Set scale at the edges of the mesh";
      return 0;
    }
  }

//  double lc_fine = lc / 2.0;
  //Waveguide to inspect:
  std::vector<int> pts;
  pts.push_back(gmsh::model::occ::addPoint(0,     0,     0, lc));    // p1
  pts.push_back(gmsh::model::occ::addPoint(L_in,  0,     0, lc));    // p2
  pts.push_back(gmsh::model::occ::addPoint(L_in,  gap_height, 0, lc));    // p3
  pts.push_back(gmsh::model::occ::addPoint(L_in + gap_width, gap_height, 0, lc));    // p4
  pts.push_back(gmsh::model::occ::addPoint(L_in + gap_width, 0,     0, lc));    // p5
  pts.push_back(gmsh::model::occ::addPoint(L_in + gap_width + iris_width, 0,     0, lc));    // p6
  pts.push_back(gmsh::model::occ::addPoint(L_in + gap_width + iris_width, gap_height, 0, lc));    // p7
  pts.push_back(gmsh::model::occ::addPoint(L_total -L_out, gap_height, 0, lc));    // p8
  pts.push_back(gmsh::model::occ::addPoint(L_total -L_out, 0,     0, lc));    // p9
  pts.push_back(gmsh::model::occ::addPoint(L_total,   0,     0, lc));    // p10
  pts.push_back(gmsh::model::occ::addPoint(L_total,   0.02,  0, lc));    // p11
  pts.push_back(gmsh::model::occ::addPoint(L_total -L_out, H,  0, lc));    // p12
  pts.push_back(gmsh::model::occ::addPoint(L_total -L_out, H-gap_height, 0, lc));    // p13
  pts.push_back(gmsh::model::occ::addPoint(L_in + gap_width + iris_width, H -gap_height, 0, lc));    // p14
  pts.push_back(gmsh::model::occ::addPoint(L_in + gap_width + iris_width, H,  0, lc));    // p15
  pts.push_back(gmsh::model::occ::addPoint(L_in + gap_width, H,  0, lc));    // p16
  pts.push_back(gmsh::model::occ::addPoint(L_in + gap_width, H-gap_height, 0, lc));    // p17
  pts.push_back(gmsh::model::occ::addPoint(L_in,   H-gap_height, 0, lc));    // p18
  pts.push_back(gmsh::model::occ::addPoint(L_in,    H, 0, lc));    // p19
  pts.push_back(gmsh::model::occ::addPoint(0,      H, 0, lc));
  // p20


  std::vector<int> lines;

  for (size_t i = 0; i < pts.size(); ++i){
    int start = pts[i];
    int end =pts[(i+1) % pts.size()];
    lines.push_back(gmsh::model::occ::addLine(start,end));
  }


    int cl = gmsh::model::occ::addCurveLoop(lines);
    int surf = gmsh::model::occ::addPlaneSurface({cl});

    gmsh::model::occ::synchronize();


    int outlet_line = -1;
    int inlet_line = -1;
    std::vector<double> line_tags_dbl;

    for(size_t i = 0; i < lines.size(); ++i) {
      double xmin, ymin, zmin, xmax, ymax, zmax;

      // Pass 6 doubles by reference
      gmsh::model::getBoundingBox(1, lines[i], xmin, ymin, zmin, xmax, ymax, zmax);

      // Check for Outlet (All x near L_total)
      if (std::abs(xmin - L_total) < 1e-6 && std::abs(xmax - L_total) < 1e-6) {
        line_tags_dbl.push_back((double)lines[i]);
          outlet_line = lines[i];
      }
      // Check for Inlet (All x near 0)
      else if (std::abs(xmin - 0.0) < 1e-6 && std::abs(xmax - 0.0) < 1e-6) {
          inlet_line = lines[i];
          line_tags_dbl.push_back((double)lines[i]);
      }



  }

// Fallback for safety
if (outlet_line == -1) {
  std::cout << "No outlet found reverting to default)" << std::endl;
    outlet_line = lines[9];
 }
if (inlet_line == -1){
  std::cout << "No inlet found reverting to default)" << std::endl;
  inlet_line = lines.back();
 }
std::vector<int> wall_lines;
  for(size_t i = 0; i < lines.size(); ++i) {
    if(lines[i] != outlet_line && lines[i] != inlet_line) {
      wall_lines.push_back(lines[i]);
      line_tags_dbl.push_back((double)lines[i]);
    }
  }

  // int f1 = gmsh::model::mesh::field::add("Distance");
  // gmsh::model::mesh::field::setNumbers(f1, "CurvesList", line_tags_dbl);
  // gmsh::model::mesh::field::setNumber(f1, "Sampling", 100);

  // int f2 = gmsh::model::mesh::field::add("Threshold");
  // gmsh::model::mesh::field::setNumber(f2,"IField", f1);

  // double size_at_wall = lc / scale;
  // double size_far_away = lc;
  // gmsh::model::mesh::field::setNumber(f2, "LcMin", size_at_wall);
  // gmsh::model::mesh::field::setNumber(f2, "LcMax", size_far_away);

  // // Distance settings
  // gmsh::model::mesh::field::setNumber(f2, "DistMin", 0.002); // Within 2mm of wall, use LcMin
  // gmsh::model::mesh::field::setNumber(f2, "DistMax", 0.01);  // After 1cm, use LcMax (interpolate in between)

  // // 3. Set this as the background mesh
  // gmsh::model::mesh::field::setAsBackgroundMesh(f2);

  // // Option: Don't extend size from boundary points (forces the Field to rule)
  // gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0);
  // =========================================================
 // gmsh::model::mesh::setRecombine(2, surf);





  gmsh::model::addPhysicalGroup(1, {inlet_line}, 1, "Inlet");
  gmsh::model::addPhysicalGroup(1, {outlet_line}, 2, "Outlet");
  gmsh::model::addPhysicalGroup(1, wall_lines, 3, "Walls");

  // Physical Surface for the domain (Attribute 1)
  gmsh::model::addPhysicalGroup(2, {surf}, 1, "Waveguide");



  gmsh::model::mesh::generate(2);
  gmsh::write(name);

  gmsh::finalize();
  return 0;
}
