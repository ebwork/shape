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

void makeRandPointsNoFaces( Eigen::MatrixXd& V, Eigen::MatrixXi& F )
{
    const int np = 30;
	// V.resize(np, 3);
    // V.Random(np, 3);
    V = Eigen::MatrixXd::Random(np, 3);

	F.resize(1, 3);
	F(0, 0) = 0;
	F(0, 1) = 1;
	F(0, 2) = 2;
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
	makeRandPointsNoFaces( V, F ); 

    // Tetrahedralize the interior
    igl::copyleft::tetgen::tetrahedralize(V,F,"", TV,TT,TF);

	// Compute barycenters
	// igl::barycenter(TV,TT,B);

    // Create faces for all of the tetrahedrons
    Eigen::MatrixXd Vtet;
    Eigen::MatrixXi Ftet;
    convertDelTets( TV, TT, Vtet, Ftet );

    // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
	viewer.data().clear();
    viewer.data().set_mesh(Vtet, Ftet);
    viewer.data().set_face_based(true);
    viewer.launch();
}
