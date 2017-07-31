package pack;

import miniufo.application.GeoFluidApplication.BoundaryCondition;
import miniufo.application.basic.VelocityFieldInSC;
import miniufo.basic.InterpolationModel.Type;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.SphericalSpatialModel;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.util.DataInterpolation;


public final class NonDivergence{
	//
	private static final String path="D:/Data/AVSIO/";
	
	//
	public static void main(String[] args){
		//interpolation();
		cStreamfunction();
	}
	
	static void cStreamfunction(){
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"VgInterp.ctl");
		DataDescriptor dd=df.getDataDescriptor();
		
		Variable[] vs=df.getVariables(new Range("",dd),"u","v");
		Variable u=vs[0];
		Variable v=vs[1];
		
		SphericalSpatialModel ssm=new SphericalSpatialModel(dd);
		VelocityFieldInSC vf=new VelocityFieldInSC(ssm);
		
		vf.setBCofX(BoundaryCondition.Periodic);
		
		Variable sf=vf.cStreamFunctionByEndlich(u,v);
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+"sf.dat");
		dw.writeData(dd,u,v,sf); dw.closeFile();
	}
	
	
	static void interpolation(){
		int denominator=9;
		
		double reso=1.0/denominator;
		
		int x=360*denominator;
		int y=180*denominator+1;
		
		System.out.println("interpolated to "+x+" * "+y+" grids with resolution of "+reso);
		
		float[] lons=new float[x];
		float[] lats=new float[y];
		
		for(int i=0;i<x;i++) lons[i]=(float)(i*reso);
		for(int j=0;j<y;j++) lats[j]=(float)(-90+j*reso);
		
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"Vg.ctl");
		
		DataInterpolation dd=new DataInterpolation(df.getDataDescriptor());
		
		dd.horizontalInterp(path+"VgInterp.dat",Type.LINEAR,lats,lons);
	}
}
