#include <iostream>
#include <Eigen/Dense>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/barycenter.h>

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
    // const int np = 30;
	// V.resize(np, 3);
    // V.Random(np, 3);
    V = Eigen::MatrixXd::Random(npoints, 3);
	for ( auto ipt=0; ipt<npoints; ++ipt)
	    for ( auto i=0; i<3; ++i)
	        V(ipt, i) = (ptmax(i) - ptmin(i)) * V(ipt, i) + ptmin(i);

	/*
	F.resize(1, 3);
	F(0, 0) = 0;
	F(0, 1) = 1;
	F(0, 2) = 2;
	*/
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
		std::vector<double> center = { -1, 1, 3} ;
		const double radius_offset = 1;
		return sqrt( (p[0] - center[0]) * (p[0] - center[0])
				   + (p[1] - center[1]) * (p[1] - center[1])
				   + (p[2] - center[2]) * (p[2] - center[2]) )
			- radius_offset;
	};
	int npoints = V.rows();
	for (unsigned ipt=0; ipt<npoints; ++ipt) {
		std::vector<double> p = { V(ipt, 0), V(ipt, 1), V(ipt, 2) };
	    sdlist.push_back( sdfunc( p) );
    }
}


int main( int argc, char *argv[] )
{
	// Tetrahedraliza set of points
	// Tetrahedralized interior
	Eigen::MatrixXd TV;
	Eigen::MatrixXi TT;
	Eigen::MatrixXi TF;
    Eigen::MatrixXd B;

	// ------------------------------------------------------
	/*
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
	// Make a simple cube
    makeCube( V, F );
    // Tetrahedralize the interior
    igl::copyleft::tetgen::tetrahedralize(V,F,"pq1.414Y", TV,TT,TF);

	// Compute barycenters
	igl::barycenter(TV,TT,B);

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
	viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    */

	// ------------------------------------------------------
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

	// Make set of points, but without any faces
	const int npoints = 100000;
	const Eigen::Vector3d shapemin(-2,  2, 4);
	const Eigen::Vector3d shapemax( 0,  4, 6);
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
    convertDelTets( TV, TTint, Vtet, Ftet );

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
	viewer.data().clear();
    viewer.data().set_mesh(Vtet, Ftet);
    viewer.data().set_face_based(true);
    viewer.launch();
}
