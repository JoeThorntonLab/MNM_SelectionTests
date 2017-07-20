/*Run YangNielsenBranchSite Null Model with delta*/

filepath = "<INPUT ALIGNMENT FILE NAME>";
outfile_path = "<OUTPUT FILE>";
LoadFunctionLibrary("../../hyphy/res/TemplateBatchFiles/TemplateModels/chooseGeneticCode.def", {"0" : "Universal"});
DataSet ds = ReadDataFile (filepath);
DataSetFilter filteredData = CreateFilter (ds,3,"","","TAA,TAG,TGA");
coding_path = LAST_FILE_PATH; /*save the last file path if you wish*/
treeString = DATAFILE_TREE;
Tree givenTree = treeString;
COUNT_GAPS_IN_FREQUENCIES = 0; /*dont want gaps/Ns to contribute to codon frequency calculations*/
HarvestFrequencies (baseFreqs,filteredData,3,1,1);

/********* Initialize variables **********/

modelKind = 0;
global P_0;
P_0:<1;
P_0:>0;
global P_1_aux;
global 	P_1 := Min(P_1_aux,1-P_0);
P_1:<1;
P_1:>0;
rateClasses = 4;
categFreqMatrix = {{P_0,P_1,(1-P_0-P_1)/(P_0+P_1)*P_0,(1-P_0-P_1)/(P_0+P_1)*P_1}} ;
categRateMatrix = {{1,2,3,4}};
category site_kind = (rateClasses, categFreqMatrix , MEAN, ,categRateMatrix, 1, 4);
global kappa_inv;
global delta;

/******** Define the evolutionary model *********/

ModelMatrixDimension = 64 - (+ _Genetic_Code["_MATRIX_ELEMENT_VALUE_==10"]);
GY_Matrix = {ModelMatrixDimension,ModelMatrixDimension};
hshift = 0;
for (h=0; h<64; h=h+1)
{
	if (_Genetic_Code[h]==10)
	{
		hshift = hshift+1;
	}
	else
	{

		vshift = hshift;
		for (v = h+1; v<64; v=v+1)
		{

                	first1  = v$16;
                	second1 = v%16$4;
                	third1  = v%4;


			first2  = h$16;
                 	second2 = h%16$4;
                 	third2  = h%4;


			diff = v-h;
			if (_Genetic_Code[v]==10)
			{
				vshift = vshift+1;
			}
			else
			{
			  	if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0) ) /* one step */
			  	{
					if (h$4==v$4)
			  		{
			  			transition = v%4;
			  			transition2= h%4;
			  		}
			  		else
			  		{
			  			if(diff%16==0)
			  			{
			  				transition = v$16;
			  				transition2= h$16;
			  			}
			  			else
			  			{
			  				transition = v%16$4;
			  				transition2= h%16$4;
			  			}
			  		}
			  		if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) /* synonymous */
			  		{
			  			if (Abs(transition-transition2)%2) /* transversion */
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := kappa_inv*synRate;
			  				GY_Matrix[v-vshift][h-hshift] := kappa_inv*synRate;
			  			}
			  			else
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := synRate;
			  				GY_Matrix[v-vshift][h-hshift] := synRate;
			  			}
				  	}
			  		else
			  		{
			  			if (Abs(transition-transition2)%2) /* transversion */
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := kappa_inv*nonSynRate;
			  				GY_Matrix[v-vshift][h-hshift] := kappa_inv*nonSynRate;
			  			}
			  			else
			  			{
			  				GY_Matrix[h-hshift][v-vshift] := nonSynRate;
			  				GY_Matrix[v-vshift][h-hshift] := nonSynRate;
			  			}
		  			}
			  	}

			else
			{

				/*populate codons with 2 differences*/

				if ( (first1!=first2) && (second1 != second2) && (third1 == third2)) 
                                {

					transition_pos1 = v$16; 
                                        transition2_pos1 = h$16;

					transition_pos2  = v%16$4; 
                                        transition2_pos2 = h%16$4;

                                        	if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) /* synonymous */
                                        	{
				                       if ( ((Abs(transition_pos1-transition2_pos1)%2) &&  (Abs(transition_pos2-transition2_pos2)%2)) || ((((transition_pos1-transition2_pos1)%2) == 0) &&  (((transition_pos2-transition2_pos2)%2) == 0)) ) 
                                                	{

								if ( (Abs(transition_pos1-transition2_pos1)%2) &&  (Abs(transition_pos2-transition2_pos2)%2) )
								{
								GY_Matrix[h-hshift][v-vshift] := synRate*delta*kappa_inv*kappa_inv;
                                                        	GY_Matrix[v-vshift][h-hshift] := synRate*delta*kappa_inv*kappa_inv;
                                                		}

                                                       		else  /* both transitions */
                                                		{
								GY_Matrix[h-hshift][v-vshift] := synRate*delta;
                                                        	GY_Matrix[v-vshift][h-hshift] := synRate*delta;

								}
							}
							else  /*transversion and transition*/
							{
   							GY_Matrix[h-hshift][v-vshift] := synRate*delta*kappa_inv;
                                                        GY_Matrix[v-vshift][h-hshift] := synRate*delta*kappa_inv;

							}

                                             }
					   else /* non-synonymous */
                                           {

    						if ( ((Abs(transition_pos1-transition2_pos1)%2) &&  (Abs(transition_pos2-transition2_pos2)%2)) || (((transition_pos1-transition2_pos1)%2 == 0) &&  ((transition_pos2-transition2_pos2)%2 == 0)) ) 
                                                {
                                                        if( (Abs(transition_pos1-transition2_pos1)%2) &&  (Abs(transition_pos2-transition2_pos2)%2))
							{
							GY_Matrix[h-hshift][v-vshift] := nonSynRate*delta*kappa_inv*kappa_inv;
                                                        GY_Matrix[v-vshift][h-hshift] := nonSynRate*delta*kappa_inv*kappa_inv;
                                                	}

                                                	else   /* both transitions */
                                                	{
                                                        GY_Matrix[h-hshift][v-vshift] := nonSynRate*delta;
                                                        GY_Matrix[v-vshift][h-hshift] := nonSynRate*delta;

                                                	}
						}

                                                else  /*transversion and transition*/
                                                {
                                                        GY_Matrix[h-hshift][v-vshift] := nonSynRate*delta*kappa_inv;
                                                        GY_Matrix[v-vshift][h-hshift] := nonSynRate*delta*kappa_inv;

                                                }

                                        }

				}

				if ( (first1 == first2) && (second1 != second2) && (third1 != third2))
                                {

							transition_pos3 = v%4; 
                                                        transition2_pos3 = h%4;

                                                        transition_pos2  = v%16$4; 
                                                        transition2_pos2 = h%16$4;


                                        if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) /* synonymous */
                                        {
                                                if ( ((Abs(transition_pos3-transition2_pos3)%2) &&  (Abs(transition_pos2-transition2_pos2)%2)) || (((transition_pos3-transition2_pos3)%2 == 0) &&  ((transition_pos2-transition2_pos2)%2 == 0)) ) 
                                                {
						 	if( (Abs(transition_pos3-transition2_pos3)%2) &&  (Abs(transition_pos2-transition2_pos2)%2) )
							{
                                                        GY_Matrix[h-hshift][v-vshift] := synRate*delta*kappa_inv*kappa_inv;
                                                        GY_Matrix[v-vshift][h-hshift] := synRate*delta*kappa_inv*kappa_inv;
                                                	}

                                                      else   /* both transitions */
                                                      {
                                                        GY_Matrix[h-hshift][v-vshift] := synRate*delta;
                                                        GY_Matrix[v-vshift][h-hshift] := synRate*delta;
						      }
						}

                                                else
                                                {
                                                        GY_Matrix[h-hshift][v-vshift] := synRate*delta*kappa_inv;
                                                        GY_Matrix[v-vshift][h-hshift] := synRate*delta*kappa_inv;

                                                }

                                        }
                                        else /* non-synonymous */
                                        {

						if ( ((Abs(transition_pos3-transition2_pos3)%2) && (Abs(transition_pos2-transition2_pos2)%2)) || (((transition_pos3-transition2_pos3)%2 == 0) &&  ((transition_pos2-transition2_pos2)%2 == 0)) )  
                                                {
							if( (Abs(transition_pos3-transition2_pos3)%2) &&  (Abs(transition_pos2-transition2_pos2)%2) )
							{
                                                        GY_Matrix[h-hshift][v-vshift] := nonSynRate*delta*kappa_inv*kappa_inv;
                                                        GY_Matrix[v-vshift][h-hshift] := nonSynRate*delta*kappa_inv*kappa_inv;
                                                	}
                                                	else /* both transitions */
                                                	{
                                                        GY_Matrix[h-hshift][v-vshift] := nonSynRate*delta;
                                                        GY_Matrix[v-vshift][h-hshift] := nonSynRate*delta;

                                                	}
						}

                                                else
                                                {
                                                        GY_Matrix[h-hshift][v-vshift] := nonSynRate*delta*kappa_inv;
                                                        GY_Matrix[v-vshift][h-hshift] := nonSynRate*delta*kappa_inv;

                                                }

                                        }
                                }

                          	if ( (first1 != first2) && (second1 == second2) && (third1 != third2)) 
                                {

    					transition_pos1 = v$16; 
                                        transition2_pos1 = h$16;

					transition_pos3 = v%4; 
                                        transition2_pos3 = h%4;


                                        if (_Genetic_Code[0][h]==_Genetic_Code[0][v]) /* synonymous */
                                        {
                                                if ( ((Abs(transition_pos1-transition2_pos1)%2) &&  (Abs(transition_pos3-transition2_pos3)%2)) || (((transition_pos1-transition2_pos1)%2 == 0) &&  ((transition_pos3-transition2_pos3)%2 == 0)) )  
						{
                                                	if( (Abs(transition_pos1-transition2_pos1)%2) &&  (Abs(transition_pos3-transition2_pos3)%2) )
							{
                                                        GY_Matrix[h-hshift][v-vshift] := synRate*delta*kappa_inv*kappa_inv;
                                                        GY_Matrix[v-vshift][h-hshift] := synRate*delta*kappa_inv*kappa_inv;
                                                	}

                                                 	else  /* both transitions */
                                                	{
                                                        	GY_Matrix[h-hshift][v-vshift] := synRate*delta;
                                                        	GY_Matrix[v-vshift][h-hshift] := synRate*delta;
                                                	}
						}

                                                else
                                                {
                                                        GY_Matrix[h-hshift][v-vshift] := synRate*delta*kappa_inv;
                                                        GY_Matrix[v-vshift][h-hshift] := synRate*delta*kappa_inv;

                                                }

                                        }
                                        else /* non-synonymous */
                                        {

                                                if ( ((Abs(transition_pos1-transition2_pos1)%2) &&  (Abs(transition_pos3-transition2_pos3)%2)) || ( ((transition_pos1-transition2_pos1)%2 == 0) &&  ( (transition_pos3-transition2_pos3)%2 == 0)) ) 
                                                {
							if((Abs(transition_pos1-transition2_pos1)%2) &&  (Abs(transition_pos3-transition2_pos3)%2) )
							{
                                                        	GY_Matrix[h-hshift][v-vshift] := nonSynRate*delta*kappa_inv*kappa_inv;
                                                        	GY_Matrix[v-vshift][h-hshift] := nonSynRate*delta*kappa_inv*kappa_inv;
                                                	}

                                                	else /* both transitions */
                                                	{
                                                        GY_Matrix[h-hshift][v-vshift] := nonSynRate*delta;
                                                        GY_Matrix[v-vshift][h-hshift] := nonSynRate*delta;

                                                	}
						}

                                                else
                                                {
                                                        GY_Matrix[h-hshift][v-vshift] := nonSynRate*delta*kappa_inv;
                                                        GY_Matrix[v-vshift][h-hshift] := nonSynRate*delta*kappa_inv;

                                                }

                                        }
                                   }

			      }

			 }
		 
		}
	}
}
PIStop = 1.0;
codonFreqs = {ModelMatrixDimension,1};
hshift = 0;
for (h=0; h<64; h=h+1)
{
	first  = h$16;
	second = h%16$4;
	third  = h%4;
	if (_Genetic_Code[h]==10)
	{
		hshift = hshift+1;
		PIStop = PIStop-baseFreqs[first][0]*baseFreqs[second][1]*baseFreqs[third][2];
		continue;
	}
	codonFreqs[h-hshift]=baseFreqs[first][0]*baseFreqs[second][1]*baseFreqs[third][2];
}
codonFreqs = codonFreqs*(1.0/PIStop);
Model GY_Model = (GY_Matrix,codonFreqs);

/*get rough starting values*/

Tree  givenTree = treeString;
LikelihoodFunction test_lf = (filteredData, givenTree);
USE_LAST_RESULTS    = 0; /*no values to re-use, start fresh for optimization*/
OPTIMIZATION_METHOD = 4; /*the only robust optimization methods are 0 and 4*/

HKY85_Matrix = {{*,t*kappa_inv,t,t*kappa_inv}
				{t*kappa_inv,*,kappa_inv*t,t}
				{t,t*kappa_inv,*,kappa_inv*t}
				{t*kappa_inv,t,kappa_inv*t,*}};


HarvestFrequencies (nucFreqs,ds,1,1,1);
Model HKY85_Model = (HKY85_Matrix,nucFreqs);
Tree nucTree = treeString;
DataSetFilter nucData = CreateFilter (ds,1);
LikelihoodFunction	nuc_lf = (nucData,nucTree);
LIKELIHOOD_FUNCTION_OUTPUT = 5;
Optimize(nuc_mle,nuc_lf);
USE_LAST_RESULTS = 1;
LIKELIHOOD_FUNCTION_OUTPUT = 3;

global omega_0;
omega_0 :< 1;

ClearConstraints (givenTree);
global omega_FG := ((site_kind==1)*omega_0+(site_kind==2)+(site_kind>2)); /* foreground model */
global omega_BG := (((site_kind==1)+(site_kind==3))*omega_0+(site_kind==2)+(site_kind==4)); /* background model */

/* apply the model to all branches */
ExecuteCommands ("givenTree."+"hg18"+".nonSynRate:=omega_FG*givenTree."+"hg18"+".synRate;");
ReplicateConstraint ("this1.?.nonSynRate:=omega_BG*this2.?.synRate",givenTree,givenTree);/*handily impose constraints for all branches*/

bNames = BranchName   (givenTree,-1);
nucBL  = BranchLength (nucTree,-1);
for (bc=0; bc<Columns(bNames)-1; bc=bc+1)
{
	ExecuteCommands ("givenTree."+bNames[bc]+".synRate=nucTree."+bNames[bc]+".t;");
}

codBL  = BranchLength (givenTree,-1);
for (bc=0; bc<Columns(bNames)-1; bc=bc+1)
{
	if (nucBL[bc]>0)
	{
		codBL[bc]=0;
		ExecuteCommands ("givenTree."+bNames[bc]+".synRate=nucTree."+bNames[bc]+".t*"+nucBL[bc]/codBL[bc]+";");
	}
}

OPTIMIZATION_PRECISION = 0.001;
LikelihoodFunction lf = (filteredData, givenTree);
while (1)
{
	Optimize 		   (mles,lf);

	LIKELIHOOD_FUNCTION_OUTPUT = 5;

	GetString 		  (lfParameters, lf, -1);
	glV 		 	= lfParameters["Local Independent"];
	stashedValues 	= {};
	for (glVI = 0; glVI < Columns (glV); glVI = glVI + 1)
	{
		ExecuteCommands ("stashedValues[\""+glV[glVI]+"\"] = " + glV[glVI] + ";\n");
	}
	glV 		 	= lfParameters["Global Independent"];
	for (glVI = 0; glVI < Columns (glV); glVI = glVI + 1)
	{
		ExecuteCommands ("stashedValues[\""+glV[glVI]+"\"] = " + glV[glVI] + ";\n");
	}

	mlBL = BranchLength (givenTree,-1);

	samples = 500;
	fprintf (outfile_path, "\nChecking for convergence by Latin Hypercube Sampling (this may take a bit of time...)\n");

	steps = 50;
	vn = {{"P_0","P_1_aux","omega_0", "omega_2"}};
	ranges	= {{0.0001,1}{0.0001,1}{0.0001,1}{1,10}};

	LFCompute (lf,LF_START_COMPUTE);

	for (sample = 0; sample < samples; sample = sample + 1)
	{
		rv = Random({1,steps}["_MATRIX_ELEMENT_COLUMN_"],0);
		for (vid = 0; vid < Columns (vn); vid = vid + 1)
		{
			ctx = vn[vid] + "=" + (ranges[vid][0] + (ranges[vid][1]-ranges[vid][0])/steps*rv[vid]);
			ExecuteCommands (ctx);
		}
		currentBL = BranchLength (givenTree,-1);
		for (bc=0; bc<Columns(bNames)-1; bc=bc+1)
		{
			if (currentBL[bc]>0)
			{
				ExecuteCommands ("givenTree."+bNames[bc]+".synRate=givenTree."+bNames[bc]+".synRate*"+mlBL[bc]/currentBL[bc]+";");
			}
		}

		LFCompute (lf,sample_value);
		if (sample_value>mles[1][0])
		{
			fprintf (outfile_path, "\nFound a better likelihood score. Restarting the optimization routine.\n");
			break;
		}
	}
	LFCompute (lf,LF_DONE_COMPUTE);

	if (sample < samples)
	{
		continue;
	}
	storedV = Rows (stashedValues);
	for (k=0; k<Columns (storedV); k=k+1)
	{
		ExecuteCommands (storedV[k] + "=" + stashedValues[storedV[k]]);
	}
	fprintf (outfile_path, "\nThe estimation procedure appears to have converged.\n");
	break;
}

LIKELIHOOD_FUNCTION_OUTPUT = 5;
fprintf (outfile_path, lf);

fprintf(outfile_path, "mles:", mles);


fprintf (outfile_path, "\nInferred rate distribution:",
		  "\n\tClass 0.  omega_0 = ", Format (omega_0, 5,3), " weight = ", Format (P_0,5,3),
		  "\n\tClass 1.  omega  := ", Format (1, 5,3), " weight = ", Format (P_1,5,3),
		  "\n\tClass 2a. Background omega_0 = ", Format (omega_0, 5,3), " foreground omega_2 = ", Format (1, 5,3), " weight = ", Format (P_0(1-P_0-P_1)/(P_0+P_1),5,3),
		  "\n\tClass 2b. Background omega  := ", Format (1, 5,3), " foreground omega_2 = ", Format (1, 5,3), " weight = ", Format (P_1(1-P_0-P_1)/(P_0+P_1),5,3), "\n"
	);
