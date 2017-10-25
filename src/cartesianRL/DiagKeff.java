//
package cartesianRL;

import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.application.contour.KeffInCTS;
import miniufo.application.statisticsModel.FilterMethods;
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
	private static final int numOfC=201;
	private static final int interpY=400;
	private static final String path=Grids.path;
	
	
	//
	public static void main(String[] args){
		String test="runH100";
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"dispInCC/"+test+"/Stat.cts");
		CtsDescriptor dd=(CtsDescriptor)(df.getDataDescriptor());
		
		Variable tr=df.getVariables(new Range("",dd),"tr6")[0];
		
		tr.replaceUndefData(-9999);
		
		ContourCartesianSpatialModel ccsm=new ContourCartesianSpatialModel(dd);
		ccsm.initContourByTracer(tr,numOfC,2,true);
		KeffInCTS keffCTS=new KeffInCTS(ccsm);
		
		Variable area   =ccsm.getAreasBoundedByContour();
		Variable tracer =ccsm.getTracerInContourCoordinate();
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
		tracer     =ccsm.interpolatedToYs(tracer     ,interpY,t);
		qGrdA      =ccsm.interpolatedToYs(qGrdA      ,interpY,t);
		intGrd2    =ccsm.interpolatedToYs(intGrd2    ,interpY,t);
		aveGrd2AlgC=ccsm.interpolatedToYs(aveGrd2AlgC,interpY,t);
		dqdye      =ccsm.interpolatedToYs(dqdye      ,interpY,t);
		Le2        =ccsm.interpolatedToYs(Le2        ,interpY,t);
		Lmin2      =ccsm.interpolatedToYs(Lmin2      ,interpY,t);
		nkeff      =ccsm.interpolatedToYs(nkeff      ,interpY,t);
		
		FilterMethods.TRunningMean(area       ,5);
		FilterMethods.TRunningMean(tracer     ,5);
		FilterMethods.TRunningMean(qGrdA      ,5);
		FilterMethods.TRunningMean(intGrd2    ,5);
		FilterMethods.TRunningMean(aveGrd2AlgC,5);
		FilterMethods.TRunningMean(dqdye      ,5);
		FilterMethods.TRunningMean(Le2        ,5);
		FilterMethods.TRunningMean(Lmin2      ,5);
		FilterMethods.TRunningMean(nkeff      ,5);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+"dispInCC/"+test+"/KeffF.dat");
		dw.writeData(dd,tracer,area,qGrdA,intGrd2,aveGrd2AlgC,dqdye,nkeff,Le2,Lmin2); dw.closeFile();
	}
}
