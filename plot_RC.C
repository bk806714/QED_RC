void read_and_draw_unp(string infile,  TGraph *gunp);
void read_and_draw_pol(string infile,  TGraph *gnoQED, TGraph *gQED);
void draw_plot_unp()
{
  // const int num_plot=3;
  const int num_plot=2;
  //const int num_plot=1;
  TGraph *gunp[10];
  string file_name1,file_name[10];
  //Int_t Q[num_plot]={3,5,10};
  Int_t  Q[num_plot]={3,10};

    for(int i=0;i<num_plot;i++)
      {
	gunp[i]= new TGraph();
	/// FOR MODIFIED STRUCTURE FUNCTIONS
	//file_name1[i]=Form("Unp_Q%d_Rot_QED_modifiedFUU",Q[i]);

	// FOR UNMODIFIED STRUCTURE FUNCTION
	//file_name1[i]=Form("Unp_Q%d_Rot_QED",Q[i]);

	// FOR SIEVERS  FUNCTION
	//file_name1[i]=Form("Siv_Q%d_Rot_QED",Q[i]);

	// for JLAB KINEMATICS
	file_name1="Unp_JLab_kine2_data";




	file_name[i]=Form("/u/home/karki/QEDSIDIS-main/%s.txt",file_name1.c_str());
	read_and_draw_unp(file_name[i], gunp[i] );
      }

 
 
  TCanvas *c2 = new TCanvas();
  c2->SetGrid();
  TLegend *legend=new TLegend(0.6,0.65,0.88,0.85);// from(x1,y1) to (x2,y2)
  legend->SetTextFont(72);
  legend->SetTextSize(0.04);
  for(int i=0;i<num_plot;i++)
    {
      gunp[i]->SetMarkerStyle(20);
      gunp[i]->GetYaxis()->SetRangeUser(0,1.8);
      gunp[i]->SetMarkerColor(i+1);
      legend->AddEntry(gunp[i],Form("Q=%d GeV",Q[i]),"P");
      // gunp3->SetLineColor(2);
      if(i==0)
	{
	  gunp[i]->Draw("AP");
	  gunp[i]->GetXaxis()->SetTitle("q_{T}/Q");
	  gunp[i]->GetYaxis()->SetTitle("#sigma_{noRC}/#sigma_{RC}");
	}
      else
	gunp[i]->Draw("Psame");
    }
  legend->Draw();
  TLine *line= new TLine(0,1,1.6,1);
  line->SetLineColor(4);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();
}







void read_and_draw_unp(string infile,  TGraph *gunp )
{

  const int num_tot=100;
  ifstream file_RC;
  file_RC.open(infile);

  if(!file_RC )
    {
    cout<<"Check your input file "<<endl;
    exit(-1);
    }
  
  Double_t qT_Q,sigma_noRC,sigma_RC,Esigma_noRC;
  Double_t RC,E_RC;
  int i=0;

  while(file_RC >> qT_Q >> sigma_noRC >>   sigma_RC >> Esigma_noRC)
 
    {
      
      RC=sigma_noRC/sigma_RC;
      //cout<<i<<"  "<<qT_Q<<"   "<<sigma_RC<<"     "<<sigma_noRC<<"    "<<RC<<endl;
      //  gunp->Set(gunp->GetN()+1);
      // gunp->SetPoint(gunp->GetN()+1,qT_Q,RC);
      gunp->SetPoint(i,qT_Q,RC);
      i++;
    
    }
   cout<<"total points "<<i<<endl;

  
   //gunp->Draw("AP");
   //gunp->SetMarkerStyle(20);
  //gunp->SetMarkerColor(2);
  // gunp10->SetLineColor(2);
  // gunp10->Draw("AP");
  //

  
}






////////////////////////For Sievers and Collins Functions //////////////////
void draw_pol( Int_t Q2=25)
{

  cout<<"Enter the Value of Q2 25 or 100"<<endl;
  cin>>Q2;
  string asymmetry="sivers";
  string kin="EIC";



  int tot=0; //only siverss
  string file_name_siv,file_name_col,file_name_tot,file_name_S,file_name_C,file_name_T;

  string file_name_mphis,file_name_pphis, file_name_m,file_name_p;
string file_name_T_p,file_name_T_m;

  //file_name_siv=Form("Siv_Q%d_Rot_QED",Q);
  file_name_siv=Form("Pol_%s_%d_Siv_constAlpha.txt",kin.c_str(),Q2);
  file_name_col=Form("Pol_%s_%d_Col_constAlpha.txt",kin.c_str(),Q2);
  //file_name_tot=Form("Pol_%s_%d.txt",Q2);
  file_name_m=Form("Pol_%s_%d_Col_Siv_mphis_constAlpha.txt",kin.c_str(),Q2);
  file_name_p=Form("Pol_%s_%d_Col_Siv_pphis_constAlpha.txt",kin.c_str(),Q2);


 file_name_S=Form("/u/home/karki/QEDSIDIS-main/newFile/%s",file_name_siv.c_str());
 file_name_C=Form("/u/home/karki/QEDSIDIS-main/newFile/%s",file_name_col.c_str());


 // file_name_T=Form("/u/home/karki/QEDSIDIS-main/newFile/%s",file_name_tot.c_str());

file_name_T_p=Form("/u/home/karki/QEDSIDIS-main/newFile/%s",file_name_p.c_str());
file_name_T_m=Form("/u/home/karki/QEDSIDIS-main/newFile/%s",file_name_m.c_str());

  ifstream file_single_asy, file_double_asy;
  int color,color_other;

  if(strcmp(asymmetry.c_str(),"sivers")==0)
   {
     file_single_asy.open(file_name_S);
     file_double_asy.open(file_name_T_m);
     color=3;
     color_other=4;
     cout<<"Analyzing Sivers Assymetry  "<<endl;
   }
  
  if(strcmp(asymmetry.c_str(),"collins")==0)
   {
     file_single_asy.open(file_name_C);
     file_double_asy.open(file_name_T_p);
     color=4;
     color_other=3;
 cout<<"Analyzing Collins Assymetry  "<<endl;
   }
  

 if(!file_single_asy)
   {
     cout<<"SINGLE ASYMMETRY FILE IS NOT AVAILABLE "<<endl;
     exit(-1);
   }
 if(!file_double_asy)
   {
     cout<<"Double asymmetry file is not available "<<endl;
     exit(-1);
   }
 const int num_tot=20;
 double qT_Q_single,sigma_noRC_single,sigma_RC_single,Esigma_noRC_single;
 double qT_Q_double,sigma_noRC_double,sigma_RC_double,Esigma_noRC_double;

 TGraph *gnoQED_single, *gQED_single;
 gnoQED_single= new TGraph();
 gQED_single= new TGraph();

 TGraph *gnoQED_double, *gQED_double;
 gnoQED_double= new TGraph();
 gQED_double= new TGraph();




 double other_asy_noRC,other_asy_RC;
 TGraph *gnoQED_other, *gQED_other;
 gnoQED_other= new TGraph();
 gQED_other= new TGraph();
for(int i=0;i<num_tot;i++)
    {
      file_single_asy >> qT_Q_single >> sigma_noRC_single >>   sigma_RC_single >> Esigma_noRC_single;
      file_double_asy >> qT_Q_double >> sigma_noRC_double >>   sigma_RC_double >> Esigma_noRC_double;

      other_asy_noRC=abs(sigma_noRC_single-sigma_noRC_double);
      other_asy_RC=abs(sigma_RC_single-sigma_RC_double);

      gnoQED_single->SetPoint(i,qT_Q_single,sigma_noRC_single);
      gQED_single->SetPoint(i,qT_Q_single,sigma_RC_single);

      gnoQED_double->SetPoint(i,qT_Q_double,sigma_noRC_double);
      gQED_double->SetPoint(i,qT_Q_double,sigma_RC_double);


      gnoQED_other->SetPoint(i,qT_Q_single,other_asy_noRC);
      gQED_other->SetPoint(i,qT_Q_single,other_asy_RC);
      // cout<<qT_Q_single<<"   "<<qT_Q_double<<"  "<<sigma_noRC_single<<"   "<<sigma_noRC_double<<"   "<<  other_asy_noRC<< "    "<<  other_asy_RC<<endl;

    }



 TCanvas *csingle;

 TLegend *legend;// from(x1,y1) to (x2,y2)

 csingle=  new TCanvas(Form("csingle%d",2), Form("csingle%d",2) );
 csingle->SetGrid();
 legend=new TLegend(0.6,0.65,0.88,0.85);
 legend->SetTextFont(72);
 legend->SetTextSize(0.04);
 
 gnoQED_single->SetLineStyle(2);
 gnoQED_single->SetLineColor(color);
 gnoQED_single->SetLineWidth(3);
 gnoQED_single->Draw("AL");


 gQED_single->Draw("Lsame");
 gQED_single->SetLineStyle(10);
 gQED_single->SetLineColor(color);
 gQED_single->SetLineWidth(3);


 

 gnoQED_other->Draw("Lsame");
 gnoQED_other->SetLineStyle(3);
 gnoQED_other->SetLineColor(color_other);
 gnoQED_other->SetLineWidth(10);



 gQED_other->Draw("Lsame");
 gQED_other->SetLineStyle(10);
 gQED_other->SetLineColor(color_other);
 gQED_other->SetLineWidth(10);


 /*
 TCanvas *cdouble;

 cdouble=  new TCanvas(Form("cdouble%d",2), Form("cdouble%d",2) );
 cdouble->SetGrid();
 gnoQED_double->SetLineStyle(3);
 gnoQED_double->SetLineColor(4);
 gnoQED_double->SetLineWidth(3);
 gnoQED_double->Draw("AL");


 gQED_double->Draw("Lsame");
 gQED_double->SetLineStyle(2);
 gQED_double->SetLineColor(4);
 gQED_double->SetLineWidth(3);
 */
}

 /*



  ifstream file_RC_siv;
  file_RC_siv.open(file_name_S);

  if(!file_RC_siv )
    {
    cout<<"Check your input file "<<endl;
    exit(-1);
    }
  
 ifstream file_RC_col;
  file_RC_col.open(file_name_C);

  if(!file_RC_col )
    {
    cout<<"Check your input COL file "<<endl;
    exit(-1);
    }
  /////////// SIVERS //////////
  Double_t qT_Q_siv,sigma_noRC_siv,sigma_RC_siv,Esigma_noRC_siv;
  TGraph *gnoQED_siv, *gQED_siv;
  gnoQED_siv= new TGraph();
  gQED_siv= new TGraph();


 
  /////////// COLLINS  //////////
  Double_t qT_Q_col,sigma_noRC_col,sigma_RC_col,Esigma_noRC_col;
  TGraph *gnoQED_col, *gQED_col;
  gnoQED_col= new TGraph();
  gQED_col= new TGraph();


  ///// TOTAL /////////////
  Double_t qT_Q_tot,sigma_noRC_tot,sigma_RC_tot,Esigma_noRC_tot;
  TGraph *gnoQED_tot, *gQED_tot;
  gnoQED_tot= new TGraph();
  gQED_tot= new TGraph();

 ifstream file_RC_tot;
  if(tot==1)
    {

  file_RC_tot.open(file_name_T);

  if(!file_RC_tot )
    {
    cout<<"Check your input TOT file "<<endl;
    exit(-1);
    }
  
    }
  // READ THE SIVERS ASYMETTRY ONLY 

  for(int i=0;i<num_tot;i++)
    {
      file_RC_siv >> qT_Q_siv >> sigma_noRC_siv >>   sigma_RC_siv >> Esigma_noRC_siv;
      file_RC_col >> qT_Q_col >> sigma_noRC_col >>   sigma_RC_col >> Esigma_noRC_col;
    





      gnoQED_siv->SetPoint(i,qT_Q_siv,sigma_noRC_siv);
      gQED_siv->SetPoint(i,qT_Q_siv,sigma_RC_siv);

      gnoQED_col->SetPoint(i,qT_Q_col,sigma_noRC_col);
      gQED_col->SetPoint(i,qT_Q_col,sigma_RC_col);
 if(tot==1)// i.e for total assymetry
   {
   file_RC_tot >> qT_Q_tot >> sigma_noRC_tot >>   sigma_RC_tot >> Esigma_noRC_tot;
      //gnoQED_tot->SetPoint(i,qT_Q_tot,sigma_noRC_tot);
      gQED_tot->SetPoint(i,qT_Q_tot,sigma_RC_tot);
   }
    }







  
  
 



  if(tot==1)
    {
  gQED_tot->Draw("Lsame");
  gQED_tot->SetLineStyle(1);
  gQED_tot->SetLineColor(2);
  gQED_tot->SetLineWidth(2);
    }
 

      // gunp3->SetLineColor(2);
      
     
//gnoQED[i]->GetYaxis()->SetRangeUser(0,1.8);
      // gnoQED->GetXaxis()->SetTitle("q_{T}/Q");
      // gnoQED->GetYaxis()->SetTitle("#sigma_{noRC}/#sigma_{RC}");
     
     
      //legend->AddEntry(gnoQED,Form("Q=%d GeV No QED ",Q),"P");
      // legend->AddEntry(gQED,Form("Q=%d GeV With QED",Q),"P");
      //legend->Draw();
    }
 
  
  TLine *line= new TLine(0,1,1.6,1);
  line->SetLineColor(4);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();
  */
  
   
  
