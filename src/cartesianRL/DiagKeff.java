//
package cartesianRL;

import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.application.contour.KeffInCTS;
import miniufo.basic.InterpolationModel.Type;
import miniufo.descriptor.CtsDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;

//
public final class DiagKeff{
	//
	private static final int numOfC=101;
	private static final int interpY=400;
	private static final String path=Grids.path;
	
	
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"Stat.cts");
		CtsDescriptor dd=(CtsDescriptor)(df.getDataDescriptor());
		
		Variable tr=df.getVariables(new Range("",dd),"tr1")[0];
		
		tr.replaceUndefData(-9999);
		
		ContourCartesianSpatialModel ccsm=new ContourCartesianSpatialModel(dd);
		ccsm.initContourByTracer(tr,10f,30f,0.1f);
		KeffInCTS keffCTS=new KeffInCTS(ccsm);
		
		Variable area   =ccsm.getAreasBoundedByContour();
		Variable grd2   =ccsm.getSquaredTracerGradient();
		Variable intGrd2=ccsm.integrateWithinContour(grd2);
		
		//ccsm.printEquivalentYs(0,0);
		
		Variable aveGrd2AlgC=keffCTS.cGradientWRTArea(intGrd2);
		Variable qGrdA      =keffCTS.cGradientWRTArea();
		Variable dqdye      =keffCTS.cDqDye();
		Variable Le2        =keffCTS.cEquivalentLengthSquare(aveGrd2AlgC,qGrdA);
		Variable Lmin2      =keffCTS.cMinimumLengthSquare();
		Variable nkeff      =keffCTS.cNormalizedKeff(aveGrd2AlgC,dqdye);
		
		/**/
		Type t=Type.LINEAR;
		
		area       =ccsm.interpolatedToYs(area       ,interpY,t);
		qGrdA      =ccsm.interpolatedToYs(qGrdA      ,interpY,t);
		intGrd2    =ccsm.interpolatedToYs(intGrd2    ,interpY,t);
		aveGrd2AlgC=ccsm.interpolatedToYs(aveGrd2AlgC,interpY,t);
		dqdye      =ccsm.interpolatedToYs(dqdye      ,interpY,t);
		Le2        =ccsm.interpolatedToYs(Le2        ,interpY,t);
		Lmin2      =ccsm.interpolatedToYs(Lmin2      ,interpY,t);
		nkeff      =ccsm.interpolatedToYs(nkeff      ,interpY,t);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+"KeffF.dat");
		dw.writeData(dd,area,qGrdA,intGrd2,aveGrd2AlgC,dqdye,nkeff,Le2,Lmin2); dw.closeFile();
	}
}
