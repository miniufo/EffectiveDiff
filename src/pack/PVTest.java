//
package pack;

import miniufo.application.contour.ContourSphericalSpatialModel;
import miniufo.application.contour.KeffInSC;
import miniufo.basic.ArrayUtil;
import miniufo.basic.InterpolationModel.Type;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import miniufo.util.DataInterpolation;


public final class PVTest{
	//
	public static void main(String[] args){
		//interpolateData(4);System.exit(0);
		
		Variable[] re1=compute(61,241);
		
		CtlDataWriteStream cdws=new CtlDataWriteStream("d:/Data/ERAInterim/Keff/PV/Keffs.dat");
		cdws.writeData(
			DiagnosisFactory.getDataDescriptor("d:/Data/ERAInterim/Keff/PV/PV.ctl"),
			ArrayUtil.concatAll(Variable.class,re1)
		); cdws.closeFile();
	}
	
	static Variable[] compute(int numOfC,int interpY){
		DiagnosisFactory df=DiagnosisFactory.parseFile("D:/Data/ERAInterim/Keff/PV/PV.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable pv=df.getVariables(new Range("",dd),"pv")[0].multiplyEq(1e6f);
		
		ContourSphericalSpatialModel cc=new ContourSphericalSpatialModel(dd);
		cc.initContourByTracer(pv,-10f,10f,0.18f);
		//cc.initContourByTracer(pv,numOfC);
		System.out.println(cc.getContours()[0][0]);
		
		KeffInSC keffSC=new KeffInSC(cc);
		
		Variable area   =cc.getAreasBoundedByContour();
		Variable intGrd2=cc.integrateWithinContour(cc.getSquaredTracerGradient());
		
		Variable aveGrd2AlgC=keffSC.cGradientWRTArea(intGrd2);
		Variable qGrdA  =keffSC.cGradientWRTArea();
		Variable dqdye  =keffSC.cDqDye();
		Variable Le2    =keffSC.cEquivalentLengthSquare(aveGrd2AlgC,qGrdA);
		Variable Lmin2  =keffSC.cMinimumLengthSquare();
		Variable nkeff  =keffSC.cNormalizedKeff(aveGrd2AlgC,dqdye);
		
		Type t=Type.LINEAR;
		
		area       =cc.interpolatedToYs(area       ,interpY,t);
		qGrdA      =cc.interpolatedToYs(qGrdA      ,interpY,t);
		intGrd2    =cc.interpolatedToYs(intGrd2    ,interpY,t);
		aveGrd2AlgC=cc.interpolatedToYs(aveGrd2AlgC,interpY,t);
		dqdye      =cc.interpolatedToYs(dqdye      ,interpY,t);
		Le2        =cc.interpolatedToYs(Le2        ,interpY,t);
		Lmin2      =cc.interpolatedToYs(Lmin2      ,interpY,t);
		nkeff      =cc.interpolatedToYs(nkeff      ,interpY,t);
		
		return new Variable[]{area,qGrdA,intGrd2,aveGrd2AlgC,dqdye,nkeff,Le2,Lmin2};
	}
	
	static void interpolateData(int res){
		DiagnosisFactory df=DiagnosisFactory.parseFile("D:/Data/ERAInterim/Keff/PV/PV.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		DataInterpolation di=new DataInterpolation(dd);
		di.horizontalInterp("d:/Data/ERAInterim/Keff/PV/PV"+res+".dat",Type.PERIODIC_CUBIC_P,Type.CUBIC_P,(dd.getYCount()-1)*res+1,dd.getXCount()*res);
	}
}
