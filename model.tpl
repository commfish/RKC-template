// --------------------------------------------------------------------------- //
//                         Size-structured model for                           //
//                               RED KING CRAB                                 //
//                        in Southeast Alaska Waters                           //
//                                                                             //
//                                VERSION 0.1                                  //
//                               November 2016                                 //
//                                                                             //
//                                  AUTHORS                                    //
//                               Katie Palof                                   //
//                         katie.palof@alaska.gov                              //
//                                                                             //
//                                                                             //
//                   Scripts in FINAL and GLOBAL sections                      //
//                         for run time statistics                             //
//                              Steven Martell                                 //
//                          martell.steve@gmail.com                            //
//                                                                             //
//                                                                             //
//                                                                             //
// --------------------------------------------------------------------------- //


DATA_SECTION

 
  // |--------------------------------------------------------------------------|
  // |MCMC OUTPUT FILE
  // |--------------------------------------------------------------------------|

     !!CLASS ofstream evalout("evalout.prj");
     !!time(&start);



  // |--------------------------------------------------------------------------|
  // | STRINGS FOR INPUT FILES                                                  
  // |--------------------------------------------------------------------------|
  // |
  // | DataFile               : data to condition the assessment model    
  // | ControlFile            : controls for years, phases, and block options 

     init_adstring DataFile;      
     init_adstring ControlFile;    

  // | BaseFileName           : file prefix used for all  model output
  // | ReportFileName         : file name to which report file is printed

     !! BaseFileName = stripExtension(DataFile);  
     !! ReportFileName = BaseFileName + adstring(".rep");

     !! cout<<"You are modeling the "<<BaseFileName<<" stock of Red King Crab"<<endl;
     !! cout<<""<<endl;



  // |--------------------------------------------------------------------------|
  // | MODEL DATA FROM CONTROL FILE                                             
  // |--------------------------------------------------------------------------|
  // | 
  // | DEBUG_FLAG             : Boolean Flag used for manual debugging

     !! ad_comm::change_datafile_name(ControlFile);

     init_int DEBUG_FLAG;

 // |-------------------------------------------------------------------------|
 // | DESIGN MATRIX FOR PARAMETER CONTROLS                                    |
 // |-------------------------------------------------------------------------|
 // | - theta_DM -> theta is a vector of estimated parameters.

    "KATIE - you'll mod this file to control whatever primary model parameters"
    "you want to model, including the parameter bounds and assumed distributions"
    "(if any)"

    init_int n_theta;
    init_matrix theta_DM(1,n_theta,1,7);
    vector    theta_ival(1,n_theta);
    vector      theta_lb(1,n_theta);
    vector      theta_ub(1,n_theta);
    ivector    theta_phz(1,n_theta);
    ivector theta_iprior(1,n_theta);
    vector      theta_p1(1,n_theta);
    vector      theta_p2(1,n_theta);
    !! theta_ival = column(theta_DM,1);
    !! theta_lb  = column(theta_DM,2);
    !! theta_ub  = column(theta_DM,3);
    !! theta_phz = ivector(column(theta_DM,4));
    !! theta_iprior = ivector(column(theta_DM,5));
    !! theta_p1 = column(theta_DM,6);
    !! theta_p2 = column(theta_DM,7);


// |---------------------------------------------------------------------------|
// | Miscellaneous Controls
// |---------------------------------------------------------------------------|
// | nMiscCont  » Number of controls to read in.
// | dMiscCont  » Vector of miscellaneous controls,


   init_int nMiscCont;
   init_vector dMiscCont(1,nMiscCont);

   number sigr;
   number sigma_sport;
   number sigma_catch;
   int     ph_rec;
   int     ph_init;
   int     ph_Fdev;
   int     ph_FdevS;
   int     ph_spr;

  // |--------------------------------------------------------------------------|
  // | END OF FILE MARKER                             
  // |--------------------------------------------------------------------------|

     init_number eof_ctl
     
     LOCAL_CALCS

     "KATIE - this simply assigns parameter names to the list of parameters on"
     "the control file"

          sigr        = dMiscCont(1);
          sigma_sport = dMiscCont(2);
          sigma_catch = dMiscCont(3);
          ph_rec      = dMiscCont(4);
          ph_init     = dMiscCont(5);
          ph_Fdev     = dMiscCont(6);
          ph_FdevS    = dMiscCont(7);
          ph_spr      = dMiscCont(8);
   

    if(eof_ctl==42) cout << BaseFileName<<".ctl has been read correctly!"<<endl;
      else 
        {    
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|      Red alert! Captain to bridge! The .ctl file is compromised!     |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof_ctl<<", but the file *should* end with 42       |"<<endl;
         cout <<"| Please check the .ctl file for errors and make sure the above calls  |"<<endl;
         cout <<"|              are matched exactly by the file's contents              |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
         exit(1); 
    }

 END_CALCS



  // |--------------------------------------------------------------------------|
  // | MODEL DATA FROM DATA FILE                                                
  // |--------------------------------------------------------------------------|
  // | This calls the data file

     !! ad_comm::change_datafile_name(DataFile);


  // |--------------------------------------------------------------------------|
  // | MODEL STRUCTURAL DIMENSIONS  - CONSTANT OVER ALL REGIONS                                            
  // |--------------------------------------------------------------------------|
  // | 
  // | nages          -> number of ages
  // | styr           -> data start year
  // | endyr          -> data end year
  // | recage         -> age at recruitment into the model
  // | spawn_fract    -> spawning month (set to May [5])
  // | srv_fract      -> survey month (set to August [8])
  // | fsh_fract       -> month when fishery takes place (set to January [1])

     init_int      nages
     init_int      styr
     init_int      endyr
     init_int      recage
     init_number   spawn_fract;
     init_number   srv_fract;
     init_number   fshy_fract;



  // |--------------------------------------------------------------------------|
  // | CALCULATED CONTAINERS                                                    
  // |--------------------------------------------------------------------------|
  // |
  // | 1) styr_rec       -> recruitment start year
  // | 2) endyr_rec      -> recruitment end year
  // |
  // | 3) endyr_sp       -> last year of estimated spawn
  // | 4) nyrs           -> number of model years
  // |
  // | 5) yy             -> vector of specific model years
  // | 6) aa             -> vector of specific model ages
  // | 7) myy            -> vector for natural mortality cubic spline
  // |
  // | 8) nrecs_est      -> number of recruitment estimates; not currently used
  // | 9) wt_mature      -> weight of mature females (for spawning biomass)
  // |                      for each region and then in toto

     int       styr_rec;
     int       styr_sp;
     
     int       endyr_sp;
     int       nyrs;
     
     ivector   yy(styr,endyr);
     ivector   aa(1,nages);
     ivector   myy(styr+1,endyr);
     
     int       nrecs_est;
     matrix    wt_mature(1,4,1,nages);




  // |--------------------------------------------------------------------------|
  // | PHYSIOLOGY  /  LIFE HISTORY                                                             
  // |--------------------------------------------------------------------------|




  // |--------------------------------------------------------------------------|
  // | TOTAL ANNUAL COMMERCIAL CATCH                                   
  // |--------------------------------------------------------------------------|


  // |--------------------------------------------------------------------------|
  // | SURVEY DATA            
  // |--------------------------------------------------------------------------|




  // |--------------------------------------------------------------------------|
  // | SIZE COMPOSITION                     
  // |--------------------------------------------------------------------------|




  // |--------------------------------------------------------------------------|
  // | CPUE                              
  // |--------------------------------------------------------------------------|




  // |--------------------------------------------------------------------------|
  // | TRANSITION MATRICES                               
  // |--------------------------------------------------------------------------|



  // |--------------------------------------------------------------------------|
  // | END OF FILE MARKER                             
  // |--------------------------------------------------------------------------|

     init_int eof_dat;
     


  // |--------------------------------------------------------------------------|
  // | MINUTIAE                              
  // |--------------------------------------------------------------------------|
  // | 
  // | 1) integer definitions 
  // | 2) 'small number' for adding to log functions
  // | 3) offsets for multinomials to allow perfect fit = 0
  
     int iyr
     int i
     int j
     int k
  
     number oo;

     number  offset;



 LOCAL_CALCS

    if(eof_dat==42) cout << BaseFileName<<".dat has been read correctly!"<<endl;
    else 
    {       
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|   ** Red alert! Captain to bridge! The .dat file is compromised! **  |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof_dat<<", but the file *should* end with 42      |"<<endl;
         cout <<"| Please check the .dat file for errors and make sure the above calls  |"<<endl;
         cout <<"|              are matched exactly by the file's contents              |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
    exit(1); 
    }


    nyrs        = endyr - styr + 1;              // Number of years for model run
    spawn_fract = (spawn_fract - 1) / 12;        // Fraction of year at which spawning occurs
    srv_fract   = (srv_fract - 1)   / 12;        // Fraction of year at which survey occurs
    fshy_fract  = (fshy_fract - 1)  / 12;        // Fraction of year at which fishery occurs
    
    wt_mature(1)   = elem_prod(wt(1),p_mature/2);      // Assumption of equal sex division
    wt_mature(2)   = elem_prod(wt(2),p_mature/2);      // Assumption of equal sex division
    wt_mature(3)   = elem_prod(wt(3),p_mature/2);      // Assumption of equal sex division
    wt_mature(4)    = elem_prod(wt(4),p_mature/2);       // Assumption of equal sex division


    yy.fill_seqadd(styr,1);
    aa.fill_seqadd(recage,1);
    myy.fill_seqadd(styr+1,1);

    oo = 0.000001;


 END_CALCS




INITIALIZATION_SECTION

	theta theta_ival;



PARAMETER_SECTION

// |---------------------------------------------------------------------------|
// | POPULATION PARAMETERS
// |---------------------------------------------------------------------------|
// | - theta(1) -> log natural mortality
// | - theta(2) -> log initial mean recruitment age-8
// | - theta(3) -> log average age-8 recruitment from styr to endyr
// | - theta(4) -> variance of log_mean_y1
// | - theta(5) -> log of mean F commercial + halibut bycatch
// | - theta(6) -> log of mean F sport fish
// | - theta(7) -> alpha (selectivity)
// | - theta(8) -> beta (selectivity)
// | - theta(9) -> delta (selectivity)
// | - theta(10) -> gamma (selectivity)
// | - theta(11) -> q catchability (commercial)
// | - theta(12) -> q catchability (iphc)

   init_bounded_number_vector theta(1,n_theta,theta_lb,theta_ub,theta_phz);
   number log_natural_mortality;
   number log_mean_rec;
   number log_mean_y1;
   number sig1;
   number log_avg_F;
   number log_avg_Fs;
   number fish_sel_slope;
   number fish_sel_a50;
   number log_q1;
   number log_q2;

   number    M;



  // |-------------------------------------------------------------------------|
  // | STOCK-RECRUITMENT                                                       |
  // |-------------------------------------------------------------------------|



  // |-------------------------------------------------------------------------|
  // | GEAR SELECTIVITY AND CATCHABILITY                                       |
  // |-------------------------------------------------------------------------|



  // |-------------------------------------------------------------------------|
  // | FISHING MORTALITY                                                       |
  // |-------------------------------------------------------------------------|



  // |-------------------------------------------------------------------------|
  // | MORTALITY                                                               |
  // |-------------------------------------------------------------------------|



  // |-------------------------------------------------------------------------|
  // | POPULATION VECTORS AND MATRICES                                         |
  // |-------------------------------------------------------------------------|



  // |-------------------------------------------------------------------------|
  // | CATCH AND SURVEY VECTORS AND MATRICES                                   |
  // |-------------------------------------------------------------------------|



  // |-------------------------------------------------------------------------|
  // | OBJECTIVE FUNCTION COMPONENTS                                           |
  // |-------------------------------------------------------------------------|


     number c_catch_like;
     number s_catch_like;
     vector surv_like(1,3);
     number age_like;

     vector penalties(1,7);
     number sprpen
     
     number Like;	
     objective_function_value obj_fun;



        
PROCEDURE_SECTION
     initializeModelParameters();
 
     Selectivity();
            if(DEBUG_FLAG == 1) cout<<"**Selectivity**"<<endl;
            
     Mortality();
            if(DEBUG_FLAG == 1) cout<<"**Mortality**"<<endl;
            
     Abundance();
            if(DEBUG_FLAG == 1) cout<<"**Abundance**"<<endl;
            
     Catch();
            if(DEBUG_FLAG == 1) cout<<"**Catch**"<<endl;
            
     Predicted();
            if(DEBUG_FLAG == 1) cout<<"**Predicted**"<<endl;
            
     Population_Summaries();
            if(DEBUG_FLAG == 1) cout<<"**Population**"<<endl; 
     
       
     Objective_Function();
     


     if (mceval_phase())
       {
         writePosteriorSamples();

         evalout<<theta<<" "
         <<log_rec_dev<<" "
         <<init_pop<<" "
         <<log_F_devs<<" "
         <<log_F_devs_sport<<" "
         <<penalties<<" "
         <<obj_fun<<" "<<endl;
       }


FUNCTION void writePosteriorSamples()
	/**
	- This function is only envoked when the -mceval
		command line option is implemented.
	*/
	//if(nf==1){
		//ofstream ofs("ssb.ps");
	//}
	//ofstream ofs0("bmT.ps",ios::app);
	//ofs0<<tot_biomassT<<endl;


       

FUNCTION void initializeModelParameters()
	//fpen = 0;

	log_natural_mortality = theta(1);
	log_mean_rec          = theta(2);
	log_mean_y1           = theta(3);
	sig1                  = theta(4);
	log_avg_F             = theta(5);
	log_avg_Fs            = theta(6);
    	fish_sel_slope        = theta(7);
    	fish_sel_a50          = theta(8);
        log_q1                = theta(9);
        log_q2                = theta(10);





FUNCTION Selectivity
  fish_sel.initialize();

  // |--------------------------------------------------------------------------|
  // | SELECTIVITY 
  // |--------------------------------------------------------------------------|

                 for (j=1;j<=nages;j++)
                   {
                     fish_sel(j)  = 1/(1+mfexp(-mfexp(fish_sel_slope)  *
                                      (j-mfexp(fish_sel_a50))));
                   }
               

                 fish_sel   =  fish_sel  / max(fish_sel);   //Scale to 1
               

FUNCTION Mortality
 // M.initialize();
  Fmort.initialize();
  Fmort2.initialize();
  F.initialize();
  F2.initialize();
  Z.initialize();
  S.initialize();


  // |--------------------------------------------------------------------------|
  // | MORTALITY 
  // |--------------------------------------------------------------------------|


  // |----------------------------|
  // | ANNUAL FISHERY DEVIATIONS
  // |----------------------------|

     Fmort   =  mfexp(log_avg_F +  log_F_devs);
            
     if(DEBUG_FLAG == 1) cout<<"Fmort"<<endl;



  // |-------------------|
  // | Fishing Mortality
  // |-------------------|

     for (j=styr; j<=endyr; j++)
       {
         F(j)  = Fmort(j)  * fish_sel;  
       }
          
     if(DEBUG_FLAG == 1) cout<<"F"<<endl;
     
      
  // |-----------------|
  // | Total Mortality
  // |-----------------|

     Z = F +  mfexp((log_natural_mortality));
     S = mfexp(-1.0*Z);

     if(DEBUG_FLAG == 1) cout<<"M"<<endl;
     
         

FUNCTION Abundance
  natageT.initialize();
  
             
     natageT(styr,1)= mfexp(log_rec_dev(styr) + log_mean_rec + (sigr*sigr)/2)/2;
     

     for(j=1;j<=nages-1;j++)
       {
         natageT(styr,j+1)=mfexp(log_mean_y1 + init_pop(j) + (sig1*sig1)/2)/2;
       }
       

     for(i=styr+1; i<=endyr; i++) // 'i' loop
       {
         natageT(i,1) = mfexp(log_rec_dev(i) + log_mean_rec + (sigr*sigr)/2)/2;            

           for(j=2; j<=nages; j++)
             {
               natageT(i,j)  = natageT(i-1,j-1)*S(i-1,j-1); 
             }                                                                         
										  
         natageT(i,nages)    += natageT(i,nages)*S(i,nages);

       } // close 'i' loop


     if(DEBUG_FLAG == 1) cout<<"Years plus"<<endl;
        


FUNCTION Catch
  catageT.initialize();
  sportcatageT.initialize();
  pred_catchT.initialize();
  pred_sportcatchT.initialize();

  // |------------------|
  // | COMMERCIAL CATCH
  // |------------------|

     for (i=styr; i<=endyr; i++)
       {
         catageT(i) = elem_div(elem_prod(elem_prod(natageT(i),F(i)),
                               (1.-S(i))),Z(i));
         pred_catchT(i) =  catageT(i)*wt(4);
       }




FUNCTION Predicted
  pred_cpueT.initialize();

  q1 = mfexp(log_q1);


  // |----------------------------------------------------------------------|
  // | COMMERCIAL FISHERY CPUE
  // |----------------------------------------------------------------------|

     for (i=1;i<=nyrs_cpueT;i++)
       {
         pred_cpueT(i) = q1* sum(elem_prod(natageT(yrs_cpueT(i)),wt(4))) / 1000; 
       }


FUNCTION Population_Summaries
  tot_biomassT.initialize();
  tot_NT.initialize();
  expl_rateT.initialize();
  spbatageT.initialize();
  spawn_biomT.initialize();
  recruitment.initialize();


     for (i=styr;i<=endyr;i++)
       {
         tot_biomassT(i) = (natageT(i)*wt(4));
         tot_NT(i)       = sum(natageT(i));
      
         expl_rateT(i)    = pred_catchT(i)/tot_biomassT(i);

         spbatageT(i)    = elem_prod(natageT(i)*mfexp(-spawn_fract*mfexp(M)),
                           wt_mature(4));
                           
         spawn_biomT(i)  = (natageT(i)*mfexp(-spawn_fract*mfexp(M))) * wt_mature(4);
         
         recruitment(i) = natageT(i,1);
      }





FUNCTION Penalties
  penalties.initialize();

    // |-----------------------------------------|
    // | Year 1 deviations penalty
    // |-----------------------------------------|
       
       penalties(1)  =norm2(init_pop+sig1*sig1/2.)/(2.*square(sig1)) +
                       size_count(init_pop)*log(sig1);



    // |-----------------------------------------|
    // | Recruitment deviations penalty
    // |-----------------------------------------|

       penalties(2)  = norm2(log_rec_dev+sigr*sigr/2.)/(2.*square(sigr)) +
                       size_count(log_rec_dev)*log(sigr);



    // |-----------------------------------------|
    // | Natural mortality penalty
    // |-----------------------------------------|
       int mpen = 0.2;
       if (last_phase())
         {
           mpen = 2;
         }
      
       penalties(3)  = 0.5*log(2*M_PI) + log(mpen) +
                        0.5*(square(log_natural_mortality - log(0.026))) / (2*square(mpen));



    // |-----------------------------------------|
    // | Stabilize F and Fs estimates
    // |-----------------------------------------|
       int fpen = 1;
       if (last_phase())
         {
           fpen = 2;
         }
      
       penalties(4)  = dnorm(log_avg_F,log(0.02),fpen);
       penalties(5)  = dnorm(log_avg_Fs,log(0.01),fpen);
       penalties(6)  = dnorm(log_mean_rec,4,fpen);
       penalties(7)  = dnorm(log_mean_y1,4,fpen);


FUNCTION Catch_Like
  c_catch_like.initialize();
  s_catch_like.initialize();

  // |----------------------------------------------------------------------|
  // | CATCH LIKELIHOODS
  // |----------------------------------------------------------------------|

     c_catch_like  +=  0.5*log(2*M_PI) + log(sigma_catch) +
                       0.5*(norm2(log(obs_catchT+oo) -
                       log(pred_catchT+oo))) / (2.*square(sigma_catch));


    s_catch_like  +=  0.5*log(2*M_PI)  +log(sigma_sport) +
                      0.5*(norm2(log(obs_sportcatchT+oo) -
                      log(pred_sportcatchT+oo))) / (2.*square(sigma_sport));



FUNCTION Surv_Like
  surv_like.initialize();
  
  // |----------------------------------------------------------------------|
  // | ROV numbers per square kilometer - log-normal
  // |----------------------------------------------------------------------|

     for (i=1; i<=nyrs_srvT; i++)
       {
         l_varT(i) = log(1. + (square(obs_srv_seT(i))/square(obs_srv_biomT(i))));


         surv_like(1) += 0.5*log(2*M_PI) + log(sqrt(l_varT(i))) +
                         0.5*(square(log(obs_srv_biomT(i)) - log(pred_srvT(i))) /
                         (2*(l_varT(i))));
       }
  // |----------------------------------------------------------------------|
  // | COMMERCIAL CPUE - normal
  // |----------------------------------------------------------------------|

     for (i=1; i<=nyrs_cpueT; i++)
       {
     
         surv_like(2) += 0.5*log(2*M_PI) + log(var_cpueT(i))+
                         0.5*(square(obs_cpueT(i)-pred_cpueT(i)))
                         / (2*var_cpueT(i)+oo);
       }


  // |----------------------------------------------------------------------|
  // | IPHC SURVEY CPUE - log-normal
  // |----------------------------------------------------------------------|

     for (i=1; i<=nyrs_cpue_iphc; i++)
       {                         

         l_var_iphcT(i) = log(1 + ((var_cpue_iphcT(i)+oo)/square(obs_cpue_iphcT(i)+oo)));

         surv_like(3)  +=  0.5*log(2*M_PI) + log(sqrt(l_var_iphcT(i))+oo) +
                           0.5*(square(log(obs_cpue_iphcT(i))-log(pred_cpue_iphcT(i))))
                           / (2*l_var_iphcT(i)+oo);
       }


FUNCTION Age_Like
  age_like.initialize();

  // |----------------------------------------------------------------------|
  // | MULTINOMIAL LIKELIHOODS FOR AGE COMPOSITIONS
  // |----------------------------------------------------------------------|

     for (i=1; i <= nyrs_fish_ageT; i++)
       {
         age_like -= nmulti_fish_ageT(i)*
                     ((oac_fishT(i) + oo) * log(age_compT(i) + oo)) ;
       }

     age_like -= offset;




FUNCTION Objective_Function
  Like.initialize();

  // |----------------------------------------------------------------------|
  // | CALL LIKELIHOOD FUNCTIONS
  // |----------------------------------------------------------------------|

     Catch_Like();
          if(DEBUG_FLAG == 1) cout<<"Catch_Like"<<endl;
     Surv_Like();
          if(DEBUG_FLAG == 1) cout<<"Surv_Like"<<endl;
     Age_Like();
          if(DEBUG_FLAG == 1) cout<<"Age_Like"<<endl;
     Penalties();
          if(DEBUG_FLAG == 1) cout<<"Penalties"<<endl;

  // |----------------------------------------------------------------------|
  // | SUM DATA LIKELIHOODS
  // |----------------------------------------------------------------------|

     Like       = c_catch_like;
     Like      += s_catch_like;
     Like      += sum(surv_like);
     Like      += age_like;

     obj_fun   = Like;         // Put here to capture the data likelihood

     obj_fun   += sum(penalties);

     



TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
  arrmblsize=39000000;

  

RUNTIME_SECTION
  maximum_function_evaluations 5000 5000  5000  5000 5000;
  convergence_criteria 0.0001;


FINAL_SECTION

   // |----------------------------------------------------------------------|
   // | Print run time statistics to the screen.
   // |----------------------------------------------------------------------|
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
        cout<<""<<endl;
	cout<<"--Objective function value: "<<obj_fun<<endl;
        cout<<""<<endl;
	cout<<"--Maximum gradient component: "<<objective_function_value::gmax<<endl;
        cout<<""<<endl;
	cout<<"*******************************************"<<endl;
	cout<<""<<endl;
 



        cout<< "O frabjous day!"<<endl;
        cout<< "The sheep frolic!"<<endl;
        cout<<""<<endl;
        cout<<"        ...moo..."<<endl;
        cout<<"            | "<<endl;
        cout<<"            | "<<endl;
        cout<<"            | "<<endl;
        cout<<"             _.%%%%%%%%%%%%%             " <<endl;
        cout<<"            //-_%%%%%%%%%%%%%            " <<endl;
        cout<<"           (_ %\\%%%%%%%%%%%%%%~             "<<endl;
        cout<<"               %%%%%%%%%%%%%%             "<<endl;
        cout<<"                 %%%%%*%%%%              "<<endl;
        cout<<"            ,,,,,,||,,,,||,,,,,          "<<endl;
        cout<<""<<endl;




GLOBALS_SECTION

	#include <admodel.h>
        #include <string.h>
        #include <string>
        #include <time.h>
        #include <sstream>
        #include <fstream>
        adstring model_name;
        adstring data_file;

	#undef REPORT
	#define REPORT(object) report << #object "\n" << setw(8) \
	<< setprecision(4) << setfixed() << object << endl;

	#undef COUT
	#define COUT(object) cout << #object "\n" << setw(6) \
	<< setprecision(3) << setfixed() << object << endl;

        time_t start,finish;
        long hour,minute,second;
        double elapsed_time;

        adstring BaseFileName;
        adstring ReportFileName;
        adstring NewFileName;

        adstring stripExtension(adstring fileName)
	 {
           /*
	     This function strips the file extension
	     from the fileName argument and returns
	     the file name without the extension.
           */

           const int length = fileName.size();
           for (int i=length; i>=0; --i)
             {
               if (fileName(i)=='.')
                 {
                   return fileName(1,i-1);
                 }
             }

           return fileName;
         }
         

 


REPORT_SECTION
 
