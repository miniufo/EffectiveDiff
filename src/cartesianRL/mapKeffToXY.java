//
package cartesianRL;

import java.util.stream.Stream;
import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.application.contour.KeffInCTS;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;

//
public final class mapKeffToXY{
	//
	private static final int numOfC=201;
	private static final DiagnosisFactory df=DiagnosisFactory.parseFile("H:/cartRL_advSchemes/Leith1_k200/Stat.cts");
	private static final DataDescriptor dd=df.getDataDescriptor();
	private static final ContourCartesianSpatialModel ccsm=new ContourCartesianSpatialModel(dd);
	
	
	//
	public static void main(String[] args){
		df.setPrinting(false);
		
		CtlDataWriteStream cdws=new CtlDataWriteStream("H:/cartRL_advSchemes/Leith1_k200/KeffXY.dat");
		cdws.setPrinting(false);
		
		df.getVariablesTimeByTime("tr1","tr2","tr3","tr4","tr5","tr6","tr7","tr8","tr9","tr10")
			.peek(vs->{
				int t=vs[0].getRange().getTRange()[0];
				if(t%100==0) System.out.println(t);
			})
			.flatMap(vs->mapAll(vs))
			.forEach(v->cdws.writeData(v));
		cdws.closeFile();
	}
	
	
	static Stream<Variable> mapAll(Variable[] vs){
		return Stream.of(vs).flatMap(v->mapToXY(v));
	}
	
	static Stream<Variable> mapToXY(Variable tr){
		ReductionDiags.changeBCToUndef(tr);
		
		ccsm.initContourByTracer(tr,numOfC,1,true);
		
		KeffInCTS keffCTS=new KeffInCTS(ccsm);
		
		Variable grd2   =ccsm.getSquaredTracerGradient();
		Variable intGrd2=ccsm.integrateWithinContour(grd2);
		
		Variable aveGrd2AlgC=keffCTS.cGradientWRTArea(intGrd2);
		Variable dqdye      =keffCTS.cDqDye();
		Variable nkeff      =keffCTS.cNormalizedKeff(aveGrd2AlgC,dqdye);
		
		Variable nkeffXY=ccsm.mapVarInContourCoordToPhysicalCoord(nkeff);
		
		return Stream.of(nkeffXY);
	}
}
