//Qengine functionalities, this also includes mkl and whatever is needed
#include <qengine/qengine.h>

//I prefer to explicitly include all the libraries I use in the header of every file ... even though qengine already includes all of them
#include <iostream>//std::cout and std::endl
#include <string>//std::string and std::string::compare
#include <chrono>//Time keeping

//Note that I have modified the cmakelist.txt to include multithreading support
#include<thread>

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
    std::cout<<"BEC optimal control comparison"<<std::endl;


    real x_range = 5;
    real x_bound = 4;
    unsigned int x_res = 512;

    //Default parameters, will be overwritten
    real g_per_N= 1.8299;
    real kinFactor = 0.36537;


    //Tweezer potential
    real sigma = 0.5;
    real amplitude = 100;
    real omega_x= 1;


    real dt = 0.0025;
    real duration = 4.0;


    std::string name = "default_time_evolve";
    //Most data is read from the coefficient json file

    std::string input;

    std::vector<count_t> Ns;

    bool GrapeBfgsL2=false;
    bool GrapeBfgsH1=false;
    bool GrapeSteepL2=false;
    bool GrapeSteepH1=false;

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
            else if ( std::string(argv[i*2+1]).compare("input")==0 )
            {
                input = std::string(argv[i*2+2]);
            }
            else if ( std::string(argv[i*2+1]).compare("output")==0 )
            {
                name = std::string(argv[i*2+2]);
            }
            else if ( std::string(argv[i*2+1]).compare("T")==0 )
            {
                real temp = stod(std::string(argv[i*2+2]));
                duration = temp;
            }
            else if ( std::string(argv[i*2+1]).compare("t_step")==0 )
            {
                real temp = stod(std::string(argv[i*2+2]));
                dt= std::abs(temp);
            }
            else if ( std::string(argv[i*2+1]).compare("use")==0 )
            {
                if ( std::string(argv[i*2+2]).compare("GrapeBFGSL2")==0 )
                    GrapeBfgsL2=true;
                else if ( std::string(argv[i*2+2]).compare("GrapeBFGSH1")==0 )
                    GrapeBfgsH1=true;
                else if ( std::string(argv[i*2+2]).compare("GrapeSteepestL2")==0 )
                    GrapeSteepL2=true;
                else if ( std::string(argv[i*2+2]).compare("GrapeSteepestH1")==0 )
                    GrapeSteepH1=true;

            }
            else if (std::string(argv[i*2+1]).compare("help")==0  )
            {
                std::cout<<"usage "<<argv[0]<<" (args)"<<std::endl;
                std::cout<<"set range to [-x_range,x_range]-------------: x_range decimal "<<std::endl;
                std::cout<<"set x-resolution----------------------------: x_res integer "<<std::endl;
                std::cout<<"set simulation time-------------------------: T decimal"<<std::endl;
                std::cout<<"set timestep--------------------------------: t_step decimal "<<std::endl;
                std::cout<<"Add N---------------------------------------: N integer "<<std::endl;
                std::cout<<"set (coefficients) input file---------------: input string"<<std::endl;
                std::cout<<"set output file-----------------------------: output string"<<std::endl;
                std::cout<<"Use optimizer-------------------------------: use Grape[BFGS/Steepest][L2/H1]"<<std::endl;
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

        g_per_N   = coefficients["g_per_N"];
        kinFactor = coefficients["kappa"];
        sigma     = coefficients["sigma"];
        amplitude = coefficients["A"];
        omega_x   = coefficients["omega_x"];

    }
    catch(...)
    {
        std::cout<<"Could not load "<<input<<std::endl;
        return 0;
    }

    std::cout<<"Working with x in [-"<<x_range<<','<<x_range<<"] with resolution "<<x_res<<std::endl;


    const auto s = gpe::makeHilbertSpace(-x_range, x_range, x_res, kinFactor);
    const auto x = s.x();


    //Set up control functions
    auto n_steps = static_cast<count_t>(duration / dt)+1;
    count_t n_print = n_steps/64;
    auto time = makeTimeControl((n_steps), dt, 0.0);

    //IF USING SHAKING INSTEAD OF TRANSFER
    //auto u = 0.1*sin(time*omega_x);//initial guess is shaking at resonance
    auto u = 2*time/duration-1;

    tell_time("Setup done ");


    DataContainer dc;
    dc["x"] =  x.vec();


    dc["use_GRAPE_BFGS_L2"]=(GrapeBfgsL2);
    dc["use_GRAPE_BFGS_H1"]=(GrapeBfgsH1);
    dc["use_GRAPE_STEEPEST_L2"]=(GrapeSteepL2);
    dc["use_GRAPE_STEEPEST_H1"]=(GrapeSteepH1);

    real sigma2 = sigma*sigma;

    auto V_dynamic=
    makePotentialFunction(
    [&x,amplitude,sigma2](const real u)
    {
        return -amplitude*(exp(-2*pow(x-u,2)/(sigma2)));
    },u.getFront().front()
    );

    //Regardless of N, the derivative of H with respect to the control is the same

    const auto dHdu = makeAnalyticDiffPotential(makePotentialFunction(
    [&x,amplitude,sigma2](const real u)
    {
        return -amplitude*4*(exp(-2*pow(x-u,2)/(sigma2)))*(x-u)/sigma2;
    },u.getFront().front()
    ));


    //Make everything which remains the same
    //Just use a very basic stopper, and a collector which only prints progress
    const auto stopper = makeStopper([](auto& optimizer) -> bool
    {
            bool stop = false;
            if (optimizer.problem().fidelity() > 0.99) {stop = true; }
            if (std::abs( optimizer.stepSize() ) < 1e-7)    {std::cout<<"Step size break"<<std::endl; stop = true;};
            if (optimizer.iteration() == 500) {stop = true;}
            return stop;
    });

    const auto maxStepSize = 5.0;
    const auto maxInitGuess = 1.0;
    const auto stepSizeFinder = makeInterpolatingStepSizeFinder(maxStepSize,maxInitGuess);


    dc["V_init"]=(V_dynamic(u.getFront()).vec());

    //Iterate over N
    count_t num = 0;


    for (count_t N : Ns)
    {
        tell_time(" N="+std::to_string(N));
        dc["N"].append(N);
        //Get wavefunctions

        auto H_dynamic = s.T()+V_dynamic+makeGpeTerm(g_per_N*N);




        tell_time("Setting up starting and target wavefunctions");
        tell_time("0/2");
        auto psi_init = makeWavefunction(H_dynamic(u.getFront())[0]);//Normalize is on by default
        tell_time("1/2");
        auto psi_trg = makeWavefunction(H_dynamic(u.getBack())[0]);
        tell_time("2/2");


        auto problem = makeStateTransferProblem(H_dynamic,dHdu,psi_init,psi_trg,u)
                     + 1e-5*Regularization(u)
                     + 1e3*Boundaries(u,RVec{-x_bound},RVec{+x_bound});




        tell_time("Setting up Grape, Steepest H1");

        count_t it_Grape_Steepest_H1 =0;
        real fid_Grape_Steepest_H1 =0;
        bool Grape_Steepest_H1_done = false;

        const auto collector_Grape_Steepest_H1 = makeCollector([&it_Grape_Steepest_H1,&fid_Grape_Steepest_H1 ](auto& optimizer) {
            it_Grape_Steepest_H1=optimizer.iteration();
            fid_Grape_Steepest_H1 = optimizer.problem().fidelity();
        });

        auto Grape_Steepest_H1 = makeGrape_bfgs_H1(problem,stopper,collector_Grape_Steepest_H1 ,stepSizeFinder);


        real Grape_Steepest_H1_time;

        //Lampda functions for running the different parameters in parrallel, alas, template lambda functions are not availible (yet) and besides, different optimizers have different requirements


        //All variables in here are modified only by one thread, so no need for any mutexes
        auto grape0 = [&GrapeSteepH1, &Grape_Steepest_H1,&collector_Grape_Steepest_H1, &Grape_Steepest_H1_time, &Grape_Steepest_H1_done ]() {
            if (GrapeSteepH1)
            {
                auto T0 = std::chrono::high_resolution_clock::now();

                collector_Grape_Steepest_H1(Grape_Steepest_H1);
                Grape_Steepest_H1.optimize();

                auto Tend = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = Tend - T0;
                Grape_Steepest_H1_time=elapsed.count();
            }
            Grape_Steepest_H1_done = true;
        };



        tell_time("Setting up Grape, Steepest L2");

        count_t it_Grape_Steepest_L2 =0;
        real fid_Grape_Steepest_L2 =0;
        bool Grape_Steepest_L2_done = false;

        const auto collector_Grape_Steepest_L2 = makeCollector([&it_Grape_Steepest_L2,&fid_Grape_Steepest_L2 ](auto& optimizer) {
            it_Grape_Steepest_L2=optimizer.iteration();
            fid_Grape_Steepest_L2 = optimizer.problem().fidelity();
        });

        auto Grape_Steepest_L2 = makeGrape_bfgs_L2(problem,stopper,collector_Grape_Steepest_L2 ,stepSizeFinder);


        real Grape_Steepest_L2_time;



        auto grape1 = [&GrapeSteepL2, &Grape_Steepest_L2,&collector_Grape_Steepest_L2, &Grape_Steepest_L2_time, &Grape_Steepest_L2_done ]() {
            if (GrapeSteepL2)
            {
                auto T0 = std::chrono::high_resolution_clock::now();

                collector_Grape_Steepest_L2(Grape_Steepest_L2);
                Grape_Steepest_L2.optimize();

                auto Tend = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = Tend - T0;
                Grape_Steepest_L2_time=elapsed.count();
            }
            Grape_Steepest_L2_done = true;
        };



        tell_time("Setting up Grape, BFGS H1");

        count_t it_Grape_BFGS_H1 =0;
        real fid_Grape_BFGS_H1 =0;
        bool Grape_BFGS_H1_done = false;

        const auto collector_Grape_BFGS_H1 = makeCollector([&it_Grape_BFGS_H1,&fid_Grape_BFGS_H1 ](auto& optimizer) {
            it_Grape_BFGS_H1=optimizer.iteration();
            fid_Grape_BFGS_H1 = optimizer.problem().fidelity();
        });

        auto Grape_BFGS_H1 = makeGrape_bfgs_H1(problem,stopper,collector_Grape_BFGS_H1 ,stepSizeFinder);


        real Grape_BFGS_H1_time;



        auto grape2 = [&GrapeBfgsH1, &Grape_BFGS_H1,&collector_Grape_BFGS_H1, &Grape_BFGS_H1_time, &Grape_BFGS_H1_done ]() {
            if (GrapeBfgsH1)
            {
                auto T0 = std::chrono::high_resolution_clock::now();

                collector_Grape_BFGS_H1(Grape_BFGS_H1);
                Grape_BFGS_H1.optimize();

                auto Tend = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = Tend - T0;
                Grape_BFGS_H1_time=elapsed.count();
            }
            Grape_BFGS_H1_done = true;
        };



        tell_time("Setting up Grape, BFGS L2");

        count_t it_Grape_BFGS_L2 =0;
        real fid_Grape_BFGS_L2 =0;
        bool Grape_BFGS_L2_done = false;

        const auto collector_Grape_BFGS_L2 = makeCollector([&it_Grape_BFGS_L2,&fid_Grape_BFGS_L2 ](auto& optimizer) {
            it_Grape_BFGS_L2=optimizer.iteration();
            fid_Grape_BFGS_L2 = optimizer.problem().fidelity();
        });

        auto Grape_BFGS_L2 = makeGrape_bfgs_L2(problem,stopper,collector_Grape_BFGS_L2 ,stepSizeFinder);


        real Grape_BFGS_L2_time;// T0grape0 = std::chrono::high_resolution_clock::now();

        //Lampda functions for running the different parameters in parrallel, alas, template lambda functions are not availible (yet) and besides, different optimizers have different requirements


        //All variables in here are modified only by one thread, so no need for any mutexes
        auto grape3 = [&GrapeBfgsL2, &Grape_BFGS_L2,&collector_Grape_BFGS_L2, &Grape_BFGS_L2_time, &Grape_BFGS_L2_done ]() {
            if (GrapeBfgsL2)
            {
                auto T0 = std::chrono::high_resolution_clock::now();

                collector_Grape_BFGS_L2(Grape_BFGS_L2);
                Grape_BFGS_L2.optimize();

                auto Tend = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = Tend - T0;
                Grape_BFGS_L2_time=elapsed.count();
            }
            Grape_BFGS_L2_done = true;
        };

        std::thread th0(grape0);
        std::thread th1(grape1);
        std::thread th2(grape2);
        std::thread th3(grape3);



        auto T0 = std::chrono::high_resolution_clock::now();
        while (!Grape_BFGS_L2_done || !Grape_BFGS_H1_done || !Grape_Steepest_L2_done || !Grape_Steepest_H1_done)
        {

            auto T1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = T1 - T0;
            if (elapsed.count()>3)//Print status every three seconds
            {
                std::cout<<std::endl;
                tell_time("----Status----");
                std::cout<<"Grape Steepest L2 "<< (Grape_Steepest_L2_done ? "not running" : "running" )<<" N= "<<it_Grape_Steepest_L2<<"/500 "<<" F="<<fid_Grape_Steepest_L2<<std::endl;
                std::cout<<"Grape Steepest H1 "<< (Grape_Steepest_H1_done ? "not running" : "running" )<<" N= "<<it_Grape_Steepest_H1<<"/500 "<<" F="<<fid_Grape_Steepest_H1<<std::endl;
                std::cout<<"Grape BFGS L2 "<< (Grape_BFGS_L2_done ? "not running" : "running" )<<" N= "<<it_Grape_BFGS_L2<<"/500 "<<" F="<<fid_Grape_BFGS_L2<<std::endl;
                std::cout<<"Grape BFGS H1 "<< (Grape_BFGS_H1_done ? "not running" : "running" )<<" N= "<<it_Grape_BFGS_H1<<"/500 "<<" F="<<fid_Grape_BFGS_H1<<std::endl;
                T0=T1;
            }
        }
        //Print status once more
        std::cout<<std::endl;
        tell_time("----Final Status----");
        std::cout<<"Grape Steepest L2 "<< (Grape_Steepest_L2_done ? "not running" : "running" )<<" N= "<<it_Grape_Steepest_L2<<"/500 "<<" F="<<fid_Grape_Steepest_L2<<std::endl;
        std::cout<<"Grape Steepest H1 "<< (Grape_Steepest_H1_done ? "not running" : "running" )<<" N= "<<it_Grape_Steepest_H1<<"/500 "<<" F="<<fid_Grape_Steepest_H1<<std::endl;
        std::cout<<"Grape BFGS L2 "<< (Grape_BFGS_L2_done ? "not running" : "running" )<<" N= "<<it_Grape_BFGS_L2<<"/500 "<<" F="<<fid_Grape_BFGS_L2<<std::endl;
        std::cout<<"Grape BFGS H1 "<< (Grape_BFGS_H1_done ? "not running" : "running" )<<" N= "<<it_Grape_BFGS_H1<<"/500 "<<" F="<<fid_Grape_BFGS_H1<<std::endl;

        //Bye bye threads
        th0.join();
        th1.join();
        th2.join();
        th3.join();



        auto u_Grape_BFGS_H1 = Grape_BFGS_H1.problem().control();
        auto u_Grape_Steepest_H1 = Grape_Steepest_H1.problem().control();
        auto u_Grape_BFGS_L2 = Grape_BFGS_L2.problem().control();
        auto u_Grape_Steepest_L2 = Grape_Steepest_L2.problem().control();

        dc["T_compute_GRAPE_BFGS_L2_"+std::to_string(num)].append(Grape_BFGS_L2_time);
        dc["Nsteps_GRAPE_BFGS_L2_"+std::to_string(num)].append(Grape_BFGS_L2.iteration());
        dc["F_GRAPE_BFGS_L2_"+std::to_string(num)].append( Grape_BFGS_L2.problem().fidelity());

        dc["T_compute_GRAPE_BFGS_H1_"+std::to_string(num)].append(Grape_BFGS_H1_time);
        dc["Nsteps_GRAPE_BFGS_H1_"+std::to_string(num)].append(Grape_BFGS_H1.iteration());
        dc["F_GRAPE_BFGS_H1_"+std::to_string(num)].append( Grape_BFGS_H1.problem().fidelity());

        dc["T_compute_GRAPE_STEEPEST_L2_"+std::to_string(num)].append(Grape_Steepest_L2_time);
        dc["Nsteps_GRAPE_STEEPEST_L2_"+std::to_string(num)].append(Grape_Steepest_L2.iteration());
        dc["F_GRAPE_STEEPEST_L2_"+std::to_string(num)].append( Grape_Steepest_L2.problem().fidelity());

        dc["T_compute_GRAPE_STEEPEST_H1_"+std::to_string(num)].append(Grape_Steepest_H1_time);
        dc["Nsteps_GRAPE_STEEPEST_H1_"+std::to_string(num)].append(Grape_Steepest_H1.iteration());
        dc["F_GRAPE_STEEPEST_H1_"+std::to_string(num)].append( Grape_Steepest_H1.problem().fidelity());

        //set up the solver to loop through the optimized versions
        auto solver0 = makeFixedTimeStepper(H_dynamic,psi_init,dt);
        auto solver1 = makeFixedTimeStepper(H_dynamic,psi_init,dt);
        auto solver2 = makeFixedTimeStepper(H_dynamic,psi_init,dt);
        auto solver3 = makeFixedTimeStepper(H_dynamic,psi_init,dt);

        dc["u_GRAPE_BFGS_L2_"+std::to_string(num)]=u_Grape_BFGS_L2.mat();
        dc["u_GRAPE_STEEPEST_L2_"+std::to_string(num)]=u_Grape_Steepest_L2.mat();
        dc["u_GRAPE_BFGS_H1_"+std::to_string(num)]=u_Grape_BFGS_H1.mat();
        dc["u_GRAPE_STEEPEST_H1_"+std::to_string(num)]=u_Grape_Steepest_H1.mat();

        double t = 0;

        //Now animate everything
        for(auto i=0; i < n_steps; i++)
        {
            //Output everything
            if (i%n_print==0)
            {
                std::cout<<"#"<<std::flush;
            }
            dc["psi2_GRAPE_BFGS_L2_"+std::to_string(num)].append(solver0.state().vec().absSquare());
            dc["psi2_GRAPE_BFGS_H1_"+std::to_string(num)].append(solver1.state().vec().absSquare());
            dc["psi2_GRAPE_STEEPEST_L2_"+std::to_string(num)].append(solver2.state().vec().absSquare());
            dc["psi2_GRAPE_STEEPEST_H1_"+std::to_string(num)].append(solver3.state().vec().absSquare());

            if (num==0)
                dc["t"].append(t);

            if(i < n_steps-1) solver0.step(u_Grape_BFGS_L2.get(i+1));
            if(i < n_steps-1) solver1.step(u_Grape_BFGS_H1.get(i+1));
            if(i < n_steps-1) solver2.step(u_Grape_Steepest_L2.get(i+1));
            if(i < n_steps-1) solver3.step(u_Grape_Steepest_H1.get(i+1));

            t += dt;
        }
        std::cout<<std::endl;
        ++num;

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
