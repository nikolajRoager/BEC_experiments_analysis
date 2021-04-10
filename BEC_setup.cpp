//Qengine functionalities, this also includes mkl and whatever is needed
#include <qengine/qengine.h>

//I prefer to explicitly include all the libraries I use in the header of every file ... even though qengine already includes all of them
#include <iostream>//std::cout and std::endl
#include <string>//std::string and std::string::compare
#include <chrono>//Time keeping

//Turns out std and qengine namespaces are in conflict, and should not be used at the same time
using namespace qengine;

std::chrono::time_point<std::chrono::high_resolution_clock> start;
//Print a message and its timestamp
void  tell_time(std::string msg)//The default value is just to ensure that the type is right
{
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    std::cout << "time = " << elapsed.count() << "s | "<<msg<<std::endl;
}

int main(int argc, char** argv)
{

    //For time keeping
    start = std::chrono::high_resolution_clock::now();
    std::cout<<"BEC setup"<<std::endl;


    //Default parameters
    real x_range = 10;
    unsigned int x_res = 512;

    real g_per_N= 1.8299;//Meanfield interaction strength
    real kinFactor = 0.36537;  // T = -kinFactor*d^2/dx^2


    //Tweezer potential
    real sigma = 0.5;
    real amplitude = 100;


    std::string name = "default";
    //Most data is read from the coefficient json file

    std::string input;

    std::vector<count_t> Ns;

    //To calculate computation time, I might want to repeat the same calculation a number of times
    count_t repeats =1;

    //I like to include the things I am editing a lot in the arguments, so that I don't have to recompile whenever I want to try something new
    for (unsigned int i = 0; i*2+1 < argc; ++i)
    {
        //Read parameters and use to set starting variables
        try
        {
            if ( std::string(argv[i*2+1]).compare("x_res")==0 )//Cast the c string to a c++ string and call the comparison method, which returns the number of differences, so 0 means a match
            {
                unsigned int temp = stoi(std::string(argv[i*2+2]));//I could also use c's atoi,
                x_res = temp;
            }
            else if ( std::string(argv[i*2+1]).compare("x_range")==0 )
            {
                real temp = stod(std::string(argv[i*2+2]));
                x_range = std::abs(temp);
            }
            else if ( std::string(argv[i*2+1]).compare("N")==0 )
            {
                auto temp = stoi(std::string(argv[i*2+2]));
                Ns.push_back(std::abs(temp));
            }
            else if ( std::string(argv[i*2+1]).compare("repeat")==0 )
            {
                auto temp = stoi(std::string(argv[i*2+2]));
                repeats = (std::abs(temp));
            }
            else if ( std::string(argv[i*2+1]).compare("input")==0 )
            {
                input = std::string(argv[i*2+2]);
            }
            else if ( std::string(argv[i*2+1]).compare("output")==0 )
            {
                name = std::string(argv[i*2+2]);
            }
            else if (std::string(argv[i*2+1]).compare("help")==0  )
            {
                std::cout<<"usage "<<argv[0]<<" (args)"<<std::endl;
                std::cout<<"set range to [-x_range,x_range]-------------: x_range integer "<<std::endl;
                std::cout<<"set x-resolution----------------------------: x_res integer "<<std::endl;
                std::cout<<"Add N---------------------------------------: N integer "<<std::endl;
                std::cout<<"set repeat-(for calculating performance)----: repeat integer "<<std::endl;
                std::cout<<"set (coefficients) input file---------------: input string"<<std::endl;
                std::cout<<"set output file-----------------------------: output string"<<std::endl;
                //End execution, so that the output doesn't get burried beneath the other output
                return 0;
            }
            else
                std::cout<<"Warning, Unknown parameter "<<argv[i*2+1]<<' '<<argv[2*i+2]<<std::endl;

        }
        catch(...)//Whatever
        {
            std::cout<<"Warning, Invalid parameter "<<argv[i*2+1]<<' '<<argv[2*i+2]<<std::endl;
        }
    }

    if (input.size()==0)
    {
        std::cout<<"Coefficients file not specified, using coefficients.json"<<std::endl;
        input = "coefficients.json";
    }
    DataContainer coefficients;

    //Default Ns
    if (Ns.size()==0)
    {
        Ns.push_back(0);//This is a trick to get g = 0; I just hope you don't try to renormalize psi later down the line
        Ns.push_back(100);
        Ns.push_back(1000);
        Ns.push_back(10000);
    }

    try
    {
        coefficients.load(input);

        g_per_N= coefficients["g_per_N"];
        kinFactor = coefficients["kappa"];
        sigma = coefficients["sigma"];
        amplitude = coefficients["A"];

    }
    catch(...)
    {
        std::cout<<"Could not load "<<input<<std::endl;
        return 0;
    }

    std::cout<<"Working with x in [-"<<x_range<<','<<x_range<<"] with resolution "<<x_res<<std::endl;


    const auto s = gpe::makeHilbertSpace(-x_range, x_range, x_res, kinFactor);
    const auto x = s.x();



    DataContainer dc;

    tell_time("Setup done");


    dc["x"] =  x.vec();


    //Tweezer potential
    auto V_static=
    makePotentialFunction(
    [&x,amplitude,sigma](real offset,real sig_scale)//Let us just use the same potential as in the rest of the experiments, while leaving this constant for now
    {
        return -amplitude*(exp(-2*pow(x-offset,2)/((sig_scale+1)*sigma*sigma)));
    },0.f,0.f);


    dc["V"]=(V_static(0,0).vec());

    //Iterate over N
    count_t num = 0;


    for (count_t i = 0; i < repeats; ++i)
        for (count_t N : Ns)
        {
            auto T0 = std::chrono::high_resolution_clock::now();
            tell_time(std::to_string(i+1)+"/"+std::to_string(repeats)+" N="+std::to_string(N));
            dc["N"].append(N);
            //Get wavefunctions

            auto H_static = s.T()+V_static+makeGpeTerm(g_per_N*N);

            auto Ground = makeWavefunction(H_static(0,0)[0]);//Normalize is on by default

            real avg_x = expectationValue(x,Ground);//Should be 0 unless something has gone wrong, but lets calculate it to be sure
            real sig2_x = expectationValue(pow(x,2),Ground)-pow(avg_x,2);
            dc["x_"+std::to_string(num)].append(avg_x);
            dc["sig2"].append(sig2_x);

            dc["psi2_"+std::to_string(num)]=(Ground.vec().absSquare());
            ++num;

            auto Tend = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> elapsed = Tend - T0;
            dc["Tcompute"].append(elapsed.count());//I want to see this
        }

    tell_time("Done, now saving ...");

    try
    {
        dc.save(name);
    }
    catch(const std::exception& E)
    {
        std::cout<<"Error saving: "<<E.what()<<std::endl;
    }

    tell_time("Done saving");

    return 0;
}
