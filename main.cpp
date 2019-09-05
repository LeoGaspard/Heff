/****************************************
 * Creator : Léo Gaspard                *
 * Mail    : leo.gaspard@outlook.fr     *
 * Git     : https://github.com/NehZio  *
 *                                      *
 * Date    : 26/08/2019                 *
 *                                      *
 * File    : main.cpp                   *
 ****************************************/


#include <iostream>
#include <iomanip>
#include <complex>
#include <Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>

int main()
{
	int norb(3);
	int nvec(3);


	Eigen::MatrixXcd c(nvec,nvec);
	Eigen::MatrixXcd d(norb,nvec);
	Eigen::MatrixXcd p(nvec,nvec);
	Eigen::MatrixXcd s(nvec,nvec);
	Eigen::MatrixXcd invSqrtS(nvec,nvec);
	Eigen::MatrixXcd pON(nvec,nvec);
	Eigen::MatrixXcd heff(nvec,nvec);
	Eigen::MatrixXd e(nvec,nvec);


	Eigen::FullPivLU<Eigen::MatrixXcd> lu;
	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eig;

	std::cout << std::fixed << std::setprecision(6);
	std::cout << std::endl << "This program will perform Heff calculations using the method described in : \n Chem. Rev. 2014, 114, 1, 429-492" << std::endl;


//*************************************************
/*
 * Manual initialization of c
 */

	c(0,1) = {-0.4867551E+00,0.0000E-10};
	c(0,0) = {-0.2653160E+00,0.0000E-10};
	c(0,2) = { 0.7520712E+00,0.0000E-10};

	c(1,1) = { 0.6126111E+00,0.0000E-10};
	c(1,0) = {-0.7456646E+00,0.0000E-10};
	c(1,2) = { 0.1330535E+00,0.0000E-10};

	c(2,1) = {-0.5748099E+00,0.0000E-10};
	c(2,0) = {-0.5748095E+00,0.0000E-10};
	c(2,2) = {-0.5748079E+00,0.0000E-10};

/*
 * Manual initializaion of e
 */
	e(0,0) =-715.10E+00;
	e(1,1) =-231.20E+00;
	e(2,2) =   0.00E+00;
//*************************************************


		

	std::cout << std::endl << "c Matrix : " << std::endl << c << std::endl;
	std::cout << std::endl << "e Matrix :  " <<  std::endl << e << std::endl << std::endl; 


	/***********************
	 * Calculation of the  *
	 * projection operator *
	 ***********************/

	std::cout << std::endl << "Projection matrix P = |I>*<I| :  " << std::endl;
        d = Eigen::MatrixXcd::Identity(norb,nvec);
	std::cout << std::endl << d << std::endl;

	/**********************
	 * Calculation of the *
	 * projected vectors  *
	 **********************/

	std::cout << std::endl << "Projected vectors into the model space :" << std::endl;
	p = c*d;
	std::cout << std::endl << p << std::endl;

	/***********************
	 * Calculation of the  *
	 * overlap matrix      *
	 ***********************/

	std::cout << std::endl << "Overlap Matrix S = c*c^T: " << std::endl;
	s = p*p.transpose();
	std::cout << std::endl << s << std::endl;

	/***********************
	 * Calculation of the  *
	 * S^(-1/2) matrix     *
	 ***********************/


	lu.compute(s.sqrt());

	if(lu.isInvertible())
	{
		invSqrtS =  lu.inverse();
	}
	else
	{
		std::cout << std::endl << "Error : The square root of the overlap matrix is not invertible" << std::endl;
		return 0;
	}

	std::cout << std::endl << "S^(-1/2) Matrix : " << std::endl;

	std::cout << std::endl << invSqrtS << std::endl;

	/***********************
	 * Calculation of the  *
	 * otrhonormalized     *
	 * projected vectors   *
	 ***********************/

	pON =  invSqrtS*p;
	std::cout << std::endl << "Orthonormalized projected vectors : " << std::endl;
	std::cout << std::endl << pON << std::endl;

	/***********************
	 * Calculation of the  *
	 * effective           *
	 * Hamiltonian         *
	 ***********************/

	heff = pON.transpose()*e*pON.conjugate();

	std::cout << std::endl << "Effective Hamiltonian |S,Ms> basis : " << std::endl;
	std::cout << std::endl << heff << std::endl;

	/***********************
	 * Diagonalization of  *
	 * the effective       *
	 * Hamiltonian         *
	 ***********************/

	eig.compute(heff);

	std::cout << std::endl << "Eigenvalues of the effective Hamiltonian :" << std::endl;
	std::cout << std::endl << eig.eigenvalues() << std::endl;

	std::cout << std::endl << "Eigenvectors of the effective Hamiltonian :" << std::endl;
	std::cout << std::endl << eig.eigenvectors() << std::endl;
	return 0;
}
