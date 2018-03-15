//
package cartesianRL;

import java.util.stream.Stream;
import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.application.contour.KeffInCTS;
import miniufo.basic.InterpolationModel.Type;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;

//
public final class DiagKeff{
	//
	private static final int numOfC=201;
	private static final int interpY=400;
	private static final String test="ctrl";
	private static final DiagnosisFactory df=DiagnosisFactory.parseFile("H:/dispInCC/"+test+"/Stat.cts");
	private static final DataDescriptor dd=df.getDataDescriptor();
	private static final ContourCartesianSpatialModel ccsm=new ContourCartesianSpatialModel(dd);
	
	
	//
	public static void main(String[] args){
		df.setPrinting(false);
		
		CtlDataWriteStream cdws=new CtlDataWriteStream("H:/dispInCC/"+test+"/test.dat"); cdws.setPrinting(false);
		df.getVariablesTimeByTime("tr1","tr2","tr3","tr4","tr5","tr6","tr7","tr8","tr9","tr10")
			.peek(vs->{
				int t=vs[0].getRange().getTRange()[0];
				if(t%100==0) System.out.println(t);
			})
			.flatMap(vs->computeKeffs(vs))
			.forEach(v->cdws.writeData(v));
		cdws.closeFile();
	}
	
	
	static Stream<Variable> computeKeffs(Variable[] vs){
		return Stream.of(vs).flatMap(v->computeKeff(v));
	}
	
	static Stream<Variable> computeKeff(Variable tr){
		ReductionDiags.changeBCToUndef(tr);
		
		ccsm.initContourByTracer(tr,numOfC,2,false,true);
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
		
		/*
		FilterMethods.TRunningMean(area       ,5);
		FilterMethods.TRunningMean(tracer     ,5);
		FilterMethods.TRunningMean(qGrdA      ,5);
		FilterMethods.TRunningMean(intGrd2    ,5);
		FilterMethods.TRunningMean(aveGrd2AlgC,5);
		FilterMethods.TRunningMean(dqdye      ,5);
		FilterMethods.TRunningMean(Le2        ,5);
		FilterMethods.TRunningMean(Lmin2      ,5);
		FilterMethods.TRunningMean(nkeff      ,5);*/
		
		return Stream.of(tracer,area,qGrdA,intGrd2,aveGrd2AlgC,dqdye,nkeff,Le2,Lmin2);
	}
}
