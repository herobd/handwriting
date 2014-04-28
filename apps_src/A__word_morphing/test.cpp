#include <stdio.h>
#include "dmorphink.h"
#include "dthresholder.h"

int main(int argc, char **argv){
	if (argc !=3)
		return 0;
	
	DMorphInk mobj;
	DImage image0;
	DImage image1;
	char stTmp[1025];
	int word0 = atoi(argv[1]);
	int word1 = atoi(argv[2]);
	
	char *stPathIn = "../../datasets/smith_jhl_vol1_lasso/smith_jhl_vol1_lasso.prj_intermediates";
	sprintf(stTmp,"%s/w_%08d.pgm",stPathIn,word0);
    	if(!image0.load(stTmp)){
      	fprintf(stderr,"couldn't load training image '%s'\n",stTmp);
      	exit(1);
    	}
    	
    	sprintf(stTmp,"%s/w_%08d.pgm",stPathIn,word1);
    	if(!image1.load(stTmp)){
      	fprintf(stderr,"couldn't load testing image '%s'\n",stTmp);
      	exit(1);
    	}
    	
    	int tval = atoi(image0.getCommentByIndex(0).c_str());
    	DThresholder::threshImage_(image0,image0, tval);
    	tval = atoi(image1.getCommentByIndex(0).c_str());
    	DThresholder::threshImage_(image1,image1, tval);

	double morphCost = mobj.getWordMorphCost(image0,
						image1,
						14,/*15 bandWidthDP*/
						0./*nonDiagonalCostDP*/,
						-1,
						-1,
						4.0,
						0.1);
	return 1;

}
