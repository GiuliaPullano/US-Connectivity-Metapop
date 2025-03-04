
//
//  simulation.cpp
//
//
//  Created by Giulia Pullano.


    
 #include "header.h"
ifstream myStream;

int main(int argc, char* argv[]){
     string output_path;
     
    
     /* GSL RANDOM NUMBER GENERATOR */
     
     unsigned long s = time(NULL)+10000;
     gsl_rng *r = gsl_rng_alloc (gsl_rng_taus);
     if (r == NULL){
         printf("allocazione con gsl fallita\n");
         exit (EXIT_FAILURE);
     }
     gsl_rng_set (r, s);
   
    //================================ DEFINE PARAMETERS =====================================//
  
    
     int n_months=12;
     
     const int n_patch= 2328;//;
     const int n_compartment=4;//

   
     const int TIME= 3000;
 
     double epsilon= 0.2703; // 3.7 days // latency
     double mu=0.3448; // 2.9 days // infecious period
     double gamma=0.2439; //(4.1 days diff) 7-days lag from cases to deaths //0.0626; //hosp+mortality rate


    //================================ DEFINE INPUT VARIABLES =====================================//
  
    
    double beta0;
    double beta=0;
    double beta1;
    
    
    string scale=argv[1];   ///spatial scale to run the model (i.e. county, state, hhs, cluster resolution)
    
    const double r0_0=stod(argv[2]);
    const double r0_1 = stod(argv[3]);

    int Nsim = 60;
    
    beta0=(double)(r0_0*mu);
    beta1=(double)(r0_1*mu);

    //================================ ALLOCATE MEMORY FOR COMPUTED VARIABLES =====================================//
    vector<int> patches_inf;
    int D_map[n_patch];
    memset(D_map, 0, sizeof D_map);
    std::vector<std::vector<int> > t0_map(n_patch, std::vector<int>(n_compartment,0));
    int S=0,E=0,I=0,R=0;
    int count=0;
    vector<vector<vector<int>>> contact_list_in(n_months, vector<vector<int>>(n_patch, vector<int>(0))); //adjacency lists incoming links
    vector<vector<vector<double>>> prob_list_in(n_months, vector<vector<double>>(n_patch, vector<double>(0))); //diffusion lists incoming links
    

    
    vector<vector<vector<int>>> contact_list_out(n_months, vector<vector<int>>(n_patch, vector<int>(0))); //adjacency lists outgoing links
    vector<vector<vector<double>>> prob_list_out(n_months, vector<vector<double>>(n_patch, vector<double>(0))); //diffusion lists outgoing links
    double prob_stay[n_months][n_patch];
    memset(prob_stay, 0, sizeof prob_stay);
    
    int  n,m;
    float prob;
    //================================ READ FILES =====================================//
  
    
    //============Read initial_condition on infected cases on March,15 =======//
    
    string initial;

    initial ="./initial_condition_1503.txt";
    ifstream read_initial(initial);
   
     
    if (read_initial.good()){
        while (read_initial >> S>> E>> I>> R) {

            t0_map[count][0] = S;
            t0_map[count][1] = E;
            t0_map[count][2] = I;
            t0_map[count][3] = R;
            count++;

        }}
    
    else cout << "No info file found 1" << endl;
    read_initial.close();

    vector<int> invasion;
    vector<int> invasion_d;


    //============Read incoming and outgoing degree in the monthly network at a given spatial scale =======//
     
     int degree_in[n_months][n_patch];
     int degree_out[n_months][n_patch];
     memset(degree_in, 0, sizeof degree_in);
     memset(degree_out, 0, sizeof degree_out);
     int month0 = 1;
     
     string  deg;
     int idx=0,in=0,out=0,month=0;
     deg="./../net_"+scale+"/degree_"+scale+".txt";
     ifstream read_degree(deg);
     if (read_degree.good()){
         while (read_degree >> idx>>month>> in>> out) {
             degree_in[month-month0][idx]=in;
             degree_out[month-month0][idx]=out;
         }}
     
     else cout << "No info file degree" << endl;
     read_degree.close();
     
     
     

   

      //==============START SIMULATIONS=============================//
    
   

     for (int isim = 0; isim < Nsim; ++isim){
       
         
         //================== ALLOCATE MEMORY FOR COMPUTED VARIABLES ============================//
        
         invasion.resize(0);
         invasion.resize(n_patch-1);
         
       
         std::vector<std::vector<int> > map_S(n_patch, std::vector<int>(TIME+1,0));
         std::vector<std::vector<int> > map_E(n_patch, std::vector<int>(TIME+1,0));
         std::vector<std::vector<int> > map_I(n_patch, std::vector<int>(TIME+1,0));
         std::vector<std::vector<int> > map_R(n_patch, std::vector<int>(TIME+1,0));
 
         std::vector<std::vector<int> > new_map_E(n_patch, std::vector<int>(TIME+1,0));
         std::vector<std::vector<int> > new_map_I(n_patch, std::vector<int>(TIME+1,0));
         std::vector<std::vector<int> > new_map_S(n_patch, std::vector<int>(TIME+1,0));
         std::vector<std::vector<int> > new_map_R(n_patch, std::vector<int>(TIME+1,0));

         std::vector<std::vector<int> > incidence(n_patch, std::vector<int>(TIME+1,0));
 
 
         patches_inf.resize(0);
         patches_inf.resize(n_patch);
         
         int month_t=0;
   
         int t0=12;
         int day=0;
         int day_tot=0;
         
         
         //================== MODEL INIZIALITATION ============================//
        
         
         for(int i = 0; i <n_patch; ++i){
             map_S[i][t0]=(int)t0_map[i][0];
             map_E[i][t0]=(int)t0_map[i][1];
             map_I[i][t0]=(int)t0_map[i][2];
             map_R[i][t0]=(int)t0_map[i][3];

             D_map[i] =map_S[i][t0]+map_E[i][t0]+map_I[i][t0]+map_R[i][t0];
             
         }
         
         for(int t= 0; t <t0; ++t){
             for(int i = 0; i <n_patch; ++i){
                 map_S[i][t]=(int)D_map[i];
                 map_E[i][t]=0;
                 map_I[i][t]=0;
                 map_R[i][t]=0;
                 map_D[i][t]=0;
                 map_D2[i][t]=0;
             }}
         


         for(int month_idx=3;month_idx<8;++month_idx){
                      
                   
                      if(month_idx==3 || month_idx==5|| month_idx==7|| month_idx==8|| month_idx==10|| month_idx==12){day_tot=31;}
                      if(month_idx==2){day_tot=28; }
                      if(month_idx==4||month_idx==11||month_idx==9){day_tot=30;}
                      if(month_idx==6){day_tot=30; }
                 
             //================== READ CONNETTIVITY NETWORK BY MONTH AND SCALE ============================//
            
                 string  matrix;
                 matrix="./../net_"+scale+"/net_monthly_matrix_sd_county_month"+ to_string(month_idx) +"_"+scale+".txt";
                 
                 cout<<"test_net: "<<month_idx<<" "<<endl;
                 //matrix="./monthly_matrix_sd_state.txt";
                 ifstream read_matrix(matrix);
                 if (read_matrix.good()){
                     while (read_matrix >> n>>m>>month>> prob) {
                         if(n!=m){
                             contact_list_out[month-month0].at(n).push_back(m);
                             prob_list_out[month-month0].at(n).push_back(prob);
                             
                             contact_list_in[month-month0].at(m).push_back(n);
                             prob_list_in[month-month0].at(m).push_back(prob);
                         }
                         if(n==m){
                             prob_stay[month-month0][n]=prob;
                             
                             
                         }}
                 }
                 
                 
                 else cout << "No info file found 2" << endl;
                 read_matrix.close();
                 
                 
             //================== RUN DAILY TIME STEPS ============================//
            

                 for(int day_idx=0;day_idx<day_tot;++day_idx){
                     day=day_idx+month_t;

                     
                     if (day<=32){
                         beta=beta0;}

                     if(day>32){
                        beta=beta1;
                         }
                     
                     
                    int degree_tot=0;
                
                     int new_time = (day + 1);
                     
                     //================== LOOP FOR EACH COUNTY ============================//
                    
                     for(int i = 0;i<n_patch; ++i){
          
                         
                         /*compute epidemic dymanic in any COUNTY */
                         
                         if(day==t0){
                             map_S[i][t0]=(int)t0_map[i][0];
                             map_E[i][t0]=(int)t0_map[i][1];
                             map_I[i][t0]=(int)t0_map[i][2];
                             map_R[i][t0]=(int)t0_map[i][3];
                             D_map[i] =map_S[i][t0]+map_E[i][t0]+map_I[i][t0]+map_R[i][t0];
                         }
                         
                         if(day<t0){
                                 map_S[i][t0]=0;
                                 map_E[i][t0]=0;
                                 map_I[i][t0]=0;
                                 map_R[i][t0]=0;
                                 map_D[i][t0]=0;
                                 map_D2[i][t0]=0;
                                 D_map[i] =map_S[i][t0]+map_E[i][t0]+map_I[i][t0]+map_R[i][t0];
                            continue;
                        
                             }
                      
                         const unint susc = map_S[i][day];
                         const unint exposed = map_E[i][day];
                         
                         const unint infect = map_I[i][day];
                         const unint recovered = map_R[i][day];
                         
          
                         const int N_sub = susc + exposed+ infect + recovered+deaths+deaths2;
                         float lambda_stay=0;
                         
                         
                         
                         
                         /* force of infection due to visitors and returing residents */
                         float  lambda_v[degree_in[month_idx-month0][i]];
                         memset(lambda_v, 0, sizeof lambda_v);
                         float  lambda_r[degree_out[month_idx-month0][i]];
                         memset(lambda_r, 0, sizeof lambda_r);
                         
                         
                         int    N_sub_j=0;
                         float  N_sub_eff=0;
                         float  N_sub_eff_j=0;
                         float  I_sub_eff_j=0;
                         float  I_sub_eff=0;
                         
                         /* compute effective popultion and infected people in any patch  */
                         N_sub_eff = prob_stay[month_idx-month0][i]*N_sub;
                         I_sub_eff = prob_stay[month_idx-month0][i]*infect;
                         
                        
                         
                         
                         int inf=infect;
                         
                         for(int j=0;j<degree_in[month_idx-month0][i];++j){
                             
                             int node =contact_list_in[month_idx-month0][i][j];
                             
                             N_sub_eff += (prob_list_in[month_idx-month0][i][j]*D_map[node]);
                             I_sub_eff += (prob_list_in[month_idx-month0][i][j]*map_I[node][day]);
                             
                         }
                         

                         //==============COMPUTE FORCE OF INFECTIN =================//
                         for(int j=0;j<degree_in[month_idx-month0][i];++j){
                             
                             int node =contact_list_in[month_idx-month0][i][j];
                             
                             lambda_v[j]=(double)((prob_stay[month_idx-month0][i]*prob_list_in[month_idx-month0][i][j]*map_I[node][day])/N_sub_eff);
                             inf+=map_I[node][day]; }
                         
                         
                         
                         int test=0;
                         for(int h=0;h<degree_out[month_idx-month0][i];++h){
                             
                             int node =contact_list_out[month_idx-month0][i][h];
                             N_sub_j=D_map[node];
                             N_sub_eff_j =prob_stay[month_idx-month0][node]*N_sub_j ;
                             I_sub_eff_j =prob_stay[month_idx-month0][node]*map_I[node][day];
                             
                             
                             
                             for(int k=0;k<degree_in[month_idx-month0][node];++k){
                                 
                                 int node2 =contact_list_in[month_idx-month0][node][k];
                                 int N_sub_k=D_map[node2];
                                 N_sub_eff_j += (prob_list_in[month_idx-month0][node][k]*N_sub_k);
                                 I_sub_eff_j += (prob_list_in[month_idx-month0][node][k]*map_I[node2][day]);
                                 inf+=map_I[node2][day];
                             }
                             lambda_r[h]=(double)((prob_list_out[month_idx-month0][i][h]*I_sub_eff_j)/N_sub_eff_j);
                             
                             test+=1;
                             
                         }
                         
    
                         
                         /* force of infection due to people not moving */
                         lambda_stay=(double)((prob_stay[month_idx-month0][i]*prob_stay[month_idx-month0][i]*infect)/N_sub_eff);
               
                         
                         
                         //==============REACTION PROCESS=================//
                         
                         int new_infected =  0;
                         int new_exposed =  0;
                         
                         int new_recovered = 0;
                         int new_deaths= 0;
                         int new_deaths2= 0;
                         degree_tot=degree_in[month_idx-month0][i]+degree_out[month_idx-month0][i];
                         
                         unint dif[degree_tot+1];
                         double diffusion_rates[degree_tot+1];
                         
                         memset(dif, 0, sizeof dif);
                         memset(diffusion_rates, 0, sizeof diffusion_rates);
                         
                         
                         int count=1;
                         double f=0;
                         
                         
                         for(int k=0;k<degree_in[month_idx-month0][i];++k){
                             
                             diffusion_rates[count]=beta*lambda_v[k];
                             
                             f+=diffusion_rates[count];
                             count+=1;
                         }
                         
                         for(int k=degree_in[month_idx-month0][i];k<degree_tot;++k){
                             int k2=k-degree_in[month_idx-month0][i];
                             diffusion_rates[count]=(beta*lambda_r[k2]);
                             f+=diffusion_rates[count];
                             count+=1;
                         }
                         
                         
                         
                         diffusion_rates[degree_tot]=(beta*lambda_stay);//last position
                         f+=diffusion_rates[degree_tot]; //last position
                    
                         diffusion_rates[0]=(1-f);
        
                         if((susc!=0) & (inf>0)){
                        
                             gsl_ran_multinomial(r,degree_tot,susc,diffusion_rates,dif);
                             
                             
                             for(int d=1;d<degree_tot+1;++d){
                              
                                 new_exposed+=dif[d]; }
               
                         }else new_exposed=0;
                         
                         
                         
                         if(exposed!=0){
                             new_infected=(int)gsl_ran_binomial(r,epsilon, exposed);
                             
                         }else new_infected=0;
                         
                   
                         if(infect!=0){
                             new_recovered=(int)gsl_ran_binomial(r,mu, infect);
                             
                         }else new_recovered=0;
                         
                            
                         
                         new_map_E[i][new_time-1] = new_exposed;
                         new_map_I[i][new_time-1] = new_infected;
                         new_map_R[i][new_time-1] = new_recovered;
                     
                         map_S[i][new_time] = (susc- new_exposed);
                         map_E[i][new_time] = (exposed+new_exposed)-new_infected;
                         map_I[i][new_time] = (infect+new_infected)-new_recovered-new_deaths;
                        
                        map_R[i][new_time] = (recovered+new_recovered);
                         
                        incidence[i][new_time] = new_infected;
                         
                         
                         if ((map_I[i][new_time]>0)| (map_E[i][new_time]>0)){
                             invasion[i]=1;
                         }
                         
                        
              
                 
            }
                       
                   }/*CLOSE DAY TIME STEP LOOP */
                 
             month_t+=day_tot;
              
                 
             
                 
             }/*CLOSE ONE SIMULATION LOOP */
         
             
        /*SAVE ONE SIMULATION OUTPUTS */
          ofstream myfile_inc;
          int n_patch2=n_patch;
          myfile_inc.open("./"+scale+"/fitted_sim_"+ to_string(isim) + "_"+ to_string(r0_3)+ "_incidence_prevalece.txt");
          if (myfile_inc.is_open()){
              for(int j=0;j<n_patch2;++j){
                  for(int t=0;t<day+2;++t){
                      myfile_inc << j <<" "<< t <<" "<<incidence[j][t]<<" "<<map_S[j][t]<<" "<<map_E[j][t]<<" "<<map_I[j][t]<<" "<<map_R[j][t]<<" "<<map_D[j][t]<<" "<<map_D2[j][t]<<" "<<D_map[j]<<"\n";
                  } }
              myfile_inc.close();
          }  else cout <<"Unable to open file3";
          
     }/*CLOSE N_SIMULATIONS LOOPS */
} /*CLOSE MAIN FUNCTION */


