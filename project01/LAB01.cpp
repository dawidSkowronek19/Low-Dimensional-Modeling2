#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <filesystem>


const double weights[8] ={
    0.1012285362903763,
    0.2223810344533745,
    0.3137066458778873,
    0.3626837833783620,
    0.3626837833783620,
    0.3137066458778873,
    0.2223810344533745,
    0.1012285362903763
};

const double integ_points[8]={
    -0.9602898564975363,
    -0.7966664774136267,
    -0.5255324099163290,
    -0.1834346424956498,
    0.1834346424956498,
    0.5255324099163290,
    0.7966664774136267,
    0.9602898564975363
};


double integrate8(double a, double b, std::function<double(double)> f)
{
    double integration=0;

    for (int alpha=1; alpha<=8; alpha++){
        double arg = (a+b)/2.0 + integ_points[alpha]*(b-a)/2.0;
        integration+=weights[alpha]*f(arg);
    }

    return integration*(b-a)/2.0;
}

int base_f_nb=0;
const int N_integ=20;
const double h = 2.0/N_integ;


double phi(double x, int i){
    return (x-1)*(x+1)*pow(x,i);
}

double rho(double x){
    return exp(-60.0*x*x);
}


double u_solution(double x, const Eigen::VectorXd& c){
    double U=0.0;
    for (int idx=0; idx<base_f_nb; idx++){
        U+=c(idx)*phi(x,idx);
    }

    return U;
}

double s_functional(const Eigen::VectorXd& c, const Eigen::VectorXd& b){
    double S=0.0;
    for (int idx=0; idx<base_f_nb; idx++){
        S+=c(idx)*b(idx);
    }

    return -0.5*S;
}

double d_functional(std::function<double(double)> u_new, std::function<double(double)> u_old){
    return sqrt(integrate8(-1.0,1.0, [u_new,u_old](double x){return (u_new(x)-u_old(x))*(u_new(x)-u_old(x));}));
}





int main(int argc, char* argv[])
{
    if (argc<2){
        return 1;
    }
    base_f_nb = std::stoi(argv[1]);
    
    Eigen::MatrixXd A(base_f_nb, base_f_nb);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(base_f_nb);



    for (int i=0; i<b.size(); i++){
        double x_start=-1.0;
        double x_end = -1.0+h;

        for (int idx=0; idx<N_integ; idx++){ //for integration
            

            b(i)+=integrate8(x_start,x_end, [i](double x){return rho(x)*phi(x,i);});
            
            x_start=x_end;
            x_end+=h;
        }

        for (int j=0; j<base_f_nb; j++){

            if ((i+j)%2==0){
                A(i,j) = 8.0*(2.0*i*j+i+j-1.0)/((i+j-1.0)*(i+j+1.0)*(i+j+3.0));
            }
            else{
                A(i,j)=0.0;
            }
        }
    }
    Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);
    //Eigen::VectorXd c = A.llt().solve(b);
    std::string filename = "u" + std::to_string(base_f_nb)+".txt";
    std::string dict = "./tmp";
    std::filesystem::create_directory(dict);

    std::ofstream u_file(dict + '/' + filename);
    
    
    double dx=0.01;
    double x_0=-1.0;
    int steps = 2.0/dx;


    for (int i=0; i<=steps; i++){
        double x_0=-1.0+i*dx;
        u_file<<x_0<<" "<<u_solution(x_0, c)<<"\n";
        x_0+=dx;
    }


    std::cout<<base_f_nb<<" | "<<u_solution(0.0, c)<<" | "<<s_functional(c, b); 
    
    

    return 0;
}