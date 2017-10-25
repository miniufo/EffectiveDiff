//
package cartesianRL;

import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.application.contour.KeffInCTS;
import miniufo.descriptor.CtsDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;

//
public final class mapKeffToXY{
	//
	private static final int numOfC=201;
	private static final String path=Grids.path;
	
	
	//
	public static void main(String[] args){
		String test="runH150";
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"dispInCC/"+test+"/Stat.cts");
		CtsDescriptor dd=(CtsDescriptor)(df.getDataDescriptor());
		ContourCartesianSpatialModel ccsm=new ContourCartesianSpatialModel(dd);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+"dispInCC/"+test+"/KeffXY.dat");
		
		for(int l=1;l<=dd.getTCount();l++){
			Variable tr=df.getVariables(new Range("t("+l+","+l+")",dd),"tr6")[0];
			
			tr.replaceUndefData(-9999);
			
			ccsm.initContourByTracer(tr,numOfC,2,true);
			
			KeffInCTS keffCTS=new KeffInCTS(ccsm);
			
			Variable grd2   =ccsm.getSquaredTracerGradient();
			Variable intGrd2=ccsm.integrateWithinContour(grd2);
			
			Variable aveGrd2AlgC=keffCTS.cGradientWRTArea(intGrd2);
			Variable dqdye      =keffCTS.cDqDye();
			Variable nkeff      =keffCTS.cNormalizedKeff(aveGrd2AlgC,dqdye);
			
			Variable nkeffXY=ccsm.mapVarInContourCoordToPhysicalCoord(nkeff);
			
			dw.writeData(nkeffXY);
		}
		
		dw.closeFile();
	}
}
