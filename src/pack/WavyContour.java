//
package pack;

import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;


public final class WavyContour{
	//
	public static void main(String[] args){
		DiagnosisFactory df=DiagnosisFactory.parseContent(
			"dset ^ModelQuater\n"+
			"title 0.25-deg resolution model\n"+
			"undef -9.99e8\n"+
			"xdef 1440 linear   0 0.25\n"+
			"ydef  361 linear   0 0.25\n"+
			"zdef    1 levels 0 1\n"+
			"tdef    1 linear 00z01Jan2000 1dy\n"+
			"vars 1\n"+
			"test 1 99 test variable\n"+
			"endvars\n"
		);
		DataDescriptor dd=df.getDataDescriptor();
		
		double resolution=1.0/4.0;
		double A=6.0,fai0=30.0,Lwn=7.0;
		
		int T=90;
		
		Variable sf=new Variable("sf",new Range(T,1,dd.getYCount(),dd.getXCount()));
		
		float[][][] sdata=sf.getData()[0];
		
		int x=sf.getXCount(),y=sf.getYCount();
		
		for(int l=0;l<T;l++)
		for(int i=0;i<x;i++){
			double lambda=i*resolution/360.0*2.0*Math.PI*Lwn+Math.PI*2.0*l/T;
			double faiThr=A*Math.sin(Math.PI*2.0*l/T)*Math.sin(lambda)+fai0;
			
			for(int j=0;j<y;j++){
				double fai=j*resolution;
				
				if(fai>faiThr) sdata[j][i][l]=1;
				else sdata[j][i][l]=0;
			}
		}
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,"d:/Data/Wavy/sf.dat");
		dw.writeData(dd,sf); dw.closeFile();
	}
}
