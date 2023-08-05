// File: genshape-del-dist/src/main.cpp
// Author: Erik Brisson
// Date created: August 2, 2023
// Purpose: Create a shape as a boundary representation as follows:
// Choose a set of points within a 3d rectangle (for now, chosen uniformaly)
// Tetrahedralize these points (note that this may add new points, if certain flags on the 
//   tetrahedalization algorithm are set
// Using a signed distance function, assign a property of "in" or "out" to each point
// Filter out the "exterior tets" (those with at least on out point)
// Create faces for the boundary of the remaining tetrahedrons
// Write an obj file using these faces

#include <iostream>
#include <sstream>

#include <Eigen/Dense>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/barycenter.h>
#include <igl/writeOBJ.h>

void makeCube( Eigen::MatrixXd& V, Eigen::MatrixXi& F )
{
  // Inline mesh of a cube
  V.resize(8,3); V <<
    0.0,0.0,0.0,
    0.0,0.0,1.0,
    0.0,1.0,0.0,
    0.0,1.0,1.0,
    1.0,0.0,0.0,
    1.0,0.0,1.0,
    1.0,1.0,0.0,
    1.0,1.0,1.0;
  F.resize(12,3); F <<
    0,6,4,
    0,2,6,
    0,3,2,
    0,1,3,
    2,7,6,
    2,3,7,
    4,6,7,
    4,7,5,
    0,4,5,
    0,5,1,
    1,5,7,
    1,7,3; // .finished();
}

void makeRandPoints( Eigen::MatrixXd& V, Eigen::MatrixXi& F,
                     const Eigen::Vector3d& ptmin, const Eigen::Vector3d& ptmax,
					 const int npoints )
{
	// Range for x, y, z = [-1, 1] in the matrix V
    V = Eigen::MatrixXd::Random(npoints, 3);
	for ( auto ipt=0; ipt<npoints; ++ipt)
	    for ( auto i=0; i<3; ++i) {
			V(ipt, i) = 0.5 * V(ipt, i) + 0.5;                       // Resize to unit cube [0,1] x [0,1] x [0,1]
	        V(ipt, i) = (ptmax(i)-ptmin(i)) * V(ipt, i) + ptmin(i);  // resize and translate to requested xample rectangle
		}
}

void convertDelTets( const Eigen::MatrixXd& TV, const Eigen::MatrixXi& TT,
                     Eigen::MatrixXd& Vtet, Eigen::MatrixXi& Ftet )
{
	int ntet = TT.rows();
	Vtet.resize(ntet*4,3);
	Ftet.resize(ntet*4,3);
	for (unsigned i=0; i<ntet;++i)
	{
		Vtet.row(i*4+0) = TV.row(TT(i,0));
		Vtet.row(i*4+1) = TV.row(TT(i,1));
		Vtet.row(i*4+2) = TV.row(TT(i,2));
		Vtet.row(i*4+3) = TV.row(TT(i,3));
		Ftet.row(i*4+0) << (i*4)+0, (i*4)+1, (i*4)+3;
		Ftet.row(i*4+1) << (i*4)+0, (i*4)+2, (i*4)+1;
		Ftet.row(i*4+2) << (i*4)+3, (i*4)+2, (i*4)+0;
		Ftet.row(i*4+3) << (i*4)+1, (i*4)+2, (i*4)+3;
	}
}

// If the value of sdlist is <= 0, it's corresponding point is on the
// interior or surface of the implicit object
// We filter out all tets whose verts are all exterior, i.e., whose sdlist > 0
void getInteriorTets( const Eigen::MatrixXi& TT, const std::vector<double>& sdlist, 
    Eigen::MatrixXi& TTint )
{
	std::vector<int> ttInteriorInd;
	int ntets = TT.rows();
	for (unsigned itet=0; itet<ntets; ++itet) {
	    bool ext = false;
		for (int i=0; i<4; ++i)
		    if (sdlist[TT(itet, i)] > 0 ) {
			    ext = true;
				break;
			}
		if (!ext)
	        ttInteriorInd.push_back(itet);
    }
	int ntetsInterior = ttInteriorInd.size();
	TTint.resize(ntetsInterior, 4);
	for (unsigned itet=0; itet<ntetsInterior; ++itet)
	    for (int i=0; i<4; ++i)
	        TTint(itet, i) = TT(ttInteriorInd[itet], i);
}

void sdist(const Eigen::MatrixXd& V, std::vector<double>& sdlist)
{
	auto sdfunc = [](std::vector<double> p) {
		std::vector<double> center = { -1, 3, 5} ;
		const double radius_offset = 1;
		return sqrt( (p[0] - center[0]) * (p[0] - center[0])
				   + (p[1] - center[1]) * (p[1] - center[1])
				   + (p[2] - center[2]) * (p[2] - center[2]) )
			- radius_offset;
	};
	int npoints = V.rows();
	for (auto ipt=0; ipt<npoints; ++ipt) {
		std::vector<double> p = { V(ipt, 0), V(ipt, 1), V(ipt, 2) };
	    sdlist.push_back( sdfunc( p) );
    }
}


void usage(const std::string& cmd, const std::string& mess)
{
    std::cout << "usage: " << cmd << "npoints fname.obj" << std::endl;
    std::cout << "exit/error message: " << mess << std::endl;
	exit(0);
}

int main( int argc, char *argv[] )
{
	if (argc != 3)
	    usage(argv[0], "too few arguments");

	std::string cmd(argv[1]);
	
	int npoints;
	std::istringstream ss(argv[1]);
	if (!(ss >> npoints)) {
	    std::cout << "Invalid number: " << argv[1] << '\n';
	    usage(cmd, "invalid \"npoints\"");
	} else if (!ss.eof()) {
	    std::cout << "Trailing characters after number: " << argv[1] << '\n';
	    usage(cmd, "invalid \"npoints\"");
	}
	
    std::string fnameOutObj(argv[2]);

	// Tetrahedraliza set of points
	// Tetrahedralized interior
	Eigen::MatrixXd TV;
	Eigen::MatrixXi TT;
	Eigen::MatrixXi TF;
    Eigen::MatrixXd B;

	// ------------------------------------------------------
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

	// Make set of points, but without any faces
	const Eigen::Vector3d shapemin(-2,  2,  4);
	const Eigen::Vector3d shapemax( 0,  4,  6);
	makeRandPoints( V, F, shapemin, shapemax, npoints ); 

    // Tetrahedralize the interior
    igl::copyleft::tetgen::tetrahedralize(V,F,"", TV,TT,TF);
    std::cout << "Number of tets in TT = " << TT.rows() << std::endl;

	// Assign distances to the points
    std::vector<double> sdlist;
    sdist(TV, sdlist);

	// Filter out the "exterior tets"
	Eigen::MatrixXi TTint;
    getInteriorTets(TT, sdlist, TTint);
    std::cout << "Number of tets in TTint = " << TTint.rows() << std::endl;

	// Compute barycenters
	// igl::barycenter(TV,TT,B);

    // Create faces for all of the tetrahedrons
    Eigen::MatrixXd Vtet;
    Eigen::MatrixXi Ftet;
    convertDelTets(TV, TTint, Vtet, Ftet);

	// Write an obj file
    igl::writeOBJ(fnameOutObj, Vtet, Ftet);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
	viewer.data().clear();
    viewer.data().set_mesh(Vtet, Ftet);
    viewer.data().set_face_based(true);
    viewer.launch();
}
