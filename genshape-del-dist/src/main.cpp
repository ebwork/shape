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
#include <fstream>
// #include <sstream>

#include <Eigen/Dense>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
// #include <igl/barycenter.h>
// #include <igl/writeOBJ.h>

void makeRandPoints( Eigen::MatrixXd& V, 
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

// Add a set of points on or outside of the bounding box,
// with signed distanc of 1
void addExteriorPoints( Eigen::MatrixXd& TV, const Eigen::Vector3d& ptmin, const Eigen::Vector3d& ptmax)
{
	int nptsAdd = 8;

    auto npts = TV.rows();
	TV.conservativeResize(npts + nptsAdd, 3);
	int k = 0;
	double eps = (1e-5)*( (ptmax[0]-ptmin[0]) + (ptmax[1]-ptmin[1]) + (ptmax[2]-ptmin[2]) );
	for (int ix=0; ix<=1; ++ix)
		for (int iy=0; iy<=1; ++iy)
			for (int iz=0; iz<=1; ++iz, ++k) {
				TV(npts+k, 0) = ix*(ptmin[0]-eps) + (1-ix)*(ptmax[0]+eps);
				TV(npts+k, 1) = iy*(ptmin[1]-eps) + (1-iy)*(ptmax[1]+eps);
				TV(npts+k, 2) = iz*(ptmin[2]-eps) + (1-iz)*(ptmax[2]+eps);
			}
}

/*
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
*/

// Just returns "interior" tetrahedrons
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
	for (auto itet=0; itet<ntetsInterior; ++itet)
	    for (auto i=0; i<4; ++i)
	        TTint(itet, i) = TT(ttInteriorInd[itet], i);
}

// Intent is to return the boundary of a set of tetrahedra
// If a tet has one exterior vert and three interior verts, we declare
// that those three verts form a boundary triangle (note that this
// may allow for the same boundary triangle to appear twice, once in
// each orientation).  Note that the %4 is needed to get boundary
// triangle correct
// If the value of sdlist is <= 0, it's corresponding point is on the
// interior or surface of the implicit object
void getBoundary( const Eigen::MatrixXi& TT, const std::vector<double>& sdlist, Eigen::MatrixXi& Ftri)
{
	std::vector<std::array<int, 3>> triBound;
	int ntets = TT.rows();
	for (auto itet=0; itet<ntets; ++itet) {
	    int next = 0;
		int iext = -1;
		for (auto i=0; i<4; ++i)
		    if (sdlist[TT(itet, i)] > 0 ) {
			    next++;
				iext = i;
			}
		if (next == 1) {
			if (iext%2)
				triBound.push_back( { TT(itet, (iext+1)%4), TT(itet, (iext+2)%4), TT(itet, (iext+3)%4) } );
			else
				triBound.push_back( { TT(itet, (iext+1)%4), TT(itet, (iext+3)%4), TT(itet, (iext+2)%4) } );
		}
	}

	auto nfaces = triBound.size();
	Ftri.resize(nfaces, 3);
	for (auto iface=0; iface<nfaces; ++iface)
	    for (auto i=0; i<3; ++i)
	        Ftri(iface, i) = triBound[iface][i];
   
}

// Convert eigen matrix of triangles to a vector of general size faces
void triToFace( const Eigen::MatrixXi& Ftri, std::vector<std::vector<int>>& faces)
{
    for (auto itri=0; itri<Ftri.rows(); ++itri)
	    faces.push_back( { Ftri(itri, 0), Ftri(itri, 1), Ftri(itri, 2) } );
}

// Signed distance function
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

// Set signed distance function to 1 for points
// outside or on the boundary of a given box
void sdistSetExterior(const Eigen::MatrixXd& V, std::vector<double>& sdlist,
	const Eigen::Vector3d& ptmin, const Eigen::Vector3d& ptmax)
{
	int npoints = V.rows();
	for (auto ipt=0; ipt<npoints; ++ipt) {
		if ( !  ( ptmin[0]<V(ipt, 0) && V(ipt, 0)<ptmin[0] )
		     && ( ptmin[1]<V(ipt, 1) && V(ipt, 1)<ptmin[1] )
		     && ( ptmin[2]<V(ipt, 2) && V(ipt, 2)<ptmin[2] ) )
	    sdlist[ipt] = 1.;
    }
}

void makeShape( const Eigen::Vector3d& shapemin, const Eigen::Vector3d& shapemax, const int npoints,
	Eigen::MatrixXd& TV, Eigen::MatrixXi& FtriBound)
{
	// Make set of points, but without any faces
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
	// makeRandPoints( V, F, shapemin, shapemax, npoints ); 
	makeRandPoints( V, shapemin, shapemax, npoints ); 

    // Add explicit exterior points
	addExteriorPoints( V, shapemin, shapemax);

    // Tetrahedralize the interior
	Eigen::MatrixXi TT;
	Eigen::MatrixXi TF;
    igl::copyleft::tetgen::tetrahedralize(V,F,"", TV,TT,TF);
    std::cout << "Number of tets in TT = " << TT.rows() << std::endl;

	// Assign distances to the points
    std::vector<double> sdlist;
    sdist(TV, sdlist);

	// Explicitly set points outside the bounding box to be exterior
    sdistSetExterior(TV, sdlist, shapemin, shapemax);

	// Get a list of triangles representing the "boundary" of the "interior" tets
	getBoundary( TT, sdlist, FtriBound);
}


// -----------------------------------------------------------------------------------
void ebWriteObj(const std::string& fname_obj, const std::string& fname_mtllib, const std::string& usemtl,
    const Eigen::MatrixXd& Vtet, const std::vector<std::vector<int>>& faces)
{
	std::ofstream myfile;
	myfile.open (fname_obj);
	if (!fname_mtllib.empty() && !usemtl.empty()) {
	    myfile << "mtllib " << fname_mtllib << std::endl;
	    myfile << "usemtl " << usemtl << std::endl;
	}
	for (auto ipt=0; ipt<Vtet.rows(); ++ipt)
	    myfile << "v " << Vtet(ipt, 0) << " " << Vtet(ipt, 1) << " " << Vtet(ipt, 2) << std::endl;
	for (auto face : faces) {
	    myfile << "f";
	    for (auto ipt : face) 
		    myfile << " " << ipt + 1;
	    myfile << std::endl;
	}
}

void ebWriteObj(const std::string& fnameBase, const int objIndex,
    const Eigen::MatrixXd& Vtet, const std::vector<std::vector<int>>& faces)
{
    ebWriteObj(fnameBase+".obj", fnameBase+".mtl", fnameBase+"_"+std::to_string(objIndex), Vtet, faces);
}

// -----------------------------------------------------------------------------------
void usage(const std::string& cmd, const std::string& mess)
{
    std::cout << "usage: " << cmd << "npoints fname_out_base" << std::endl;
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
	
    std::string fnameOutBase(argv[2]);

	Eigen::MatrixXd TV;
	Eigen::MatrixXi FtriBound;

	// Using the bounds specified by shapemin and shapmax, choose npoints within
	// the bounds, then use a signed distance function to iedenify "interior"
	// and "exterior" points.  Use tetrahalization to create a polyhedral shape.
	const Eigen::Vector3d shapemin(-2,  2,  4);
	const Eigen::Vector3d shapemax( 0,  4,  6);
    makeShape( shapemin, shapemax, npoints, TV, FtriBound);

    // Convert triangle matrix to general face form
	std::vector<std::vector<int>> boundary;
	triToFace(FtriBound, boundary);

	std::cout << "FtriBound.rows() = " << FtriBound.rows() << std::endl;
	std::cout << "boundary.size() = " << boundary.size() << std::endl;

	// Write an obj file using my own routine
    ebWriteObj(fnameOutBase, 2, TV, boundary);

	// Write an obj file
    // igl::writeOBJ(fnameOutObj, Vtet, Ftet);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
	viewer.data().clear();
    // viewer.data().set_mesh(Vtet, Ftet);
    viewer.data().set_mesh(TV, FtriBound);
    viewer.data().set_face_based(true);
    viewer.launch();
}
