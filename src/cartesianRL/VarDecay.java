//
package cartesianRL;

import java.util.concurrent.TimeUnit;
import miniufo.application.contour.ContourCartesianSpatialModel;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;
import miniufo.util.TicToc;


//
public final class VarDecay{
	//
	private static final int y=Parameters.y;
	private static final int x=Parameters.x;
	
	private static final String path=Parameters.path;
	
	//
	public static void main(String[] args){
		int trNum=16;
		
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"Stat.cts");
		DataDescriptor dd=df.getDataDescriptor();
		
		Range r=new Range("",dd);
		
		Variable[] vars =new Variable[trNum*4];
		Variable[] areas=new Variable[trNum];
		
		for(int m=0,mm=0;m<trNum;m++){
			TicToc.tic("start computing tr"+(m+1));
			float mean=2+m; if(mean>9) mean-=8;
			
			Variable tr=df.getVariables(r,"tr"+(m+1))[0];
			
			Variable[] tmp=cMeanVarianceAndExtremes(tr);
			vars[mm++]=tmp[0];
			vars[mm++]=tmp[1];
			vars[mm++]=tmp[2];
			vars[mm++]=tmp[3];
			
			areas[m]=cContourArea(tr,dd,mean+1,mean-1);
			areas[m].setName("area"+(m+1));
			TicToc.toc(TimeUnit.SECONDS);
		}
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+"VarDecayF.dat");
		dw.writeData(dd,vars); dw.closeFile();
		
		dw=DataIOFactory.getDataWrite(dd,path+"AreasF.dat");
		dw.writeData(dd,areas); dw.closeFile();
	}
	
	static Variable cContourArea(Variable tr,DataDescriptor dd,float csouth,float cnorth){
		tr.replaceUndefData(-9999);
		
		ContourCartesianSpatialModel ccsm=new ContourCartesianSpatialModel(dd);
		ccsm.initContourByTracer(tr,csouth,cnorth,-0.25f,1,false);
		
		Variable area=ccsm.getAreasBoundedByContour();
		
		//Type t=Type.LINEAR;
		
		//area=ccsm.interpolatedToYs(area,100,t);
		
		  tr.replaceUndefData(0);
		area.replaceUndefData(0);
		
		return area;
	}
	
	static Variable[] cMeanVarianceAndExtremes(Variable v){
		int t=v.getTCount();
		float undef=v.getUndef();
		
		Variable ave=new Variable(v.getName()+"ave",new Range(t,1,1,1));
		Variable var=new Variable(v.getName()+"var",new Range(t,1,1,1));
		Variable max=new Variable(v.getName()+"max",new Range(t,1,1,1));
		Variable min=new Variable(v.getName()+"min",new Range(t,1,1,1));
		
		ave.setCommentAndUnit("average  of "+v.getName()+" (unit)"  );
		var.setCommentAndUnit("variance of "+v.getName()+" (unit^2)");
		max.setCommentAndUnit("maximum  of "+v.getName()+" (unit)"  );
		min.setCommentAndUnit("minimum  of "+v.getName()+" (unit)"  );
		
		float[][][] vdata=v.getData()[0];
		float[] avedata=ave.getData()[0][0][0];
		float[] vardata=var.getData()[0][0][0];
		float[] maxdata=max.getData()[0][0][0];
		float[] mindata=min.getData()[0][0][0];
		
		for(int l=0;l<t;l++){
			float maxV=Float.MIN_VALUE;
			float minV=Float.MAX_VALUE;
			
			double sumVar=0,sumAve=0; int count=0;
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(vdata[j][i][l]!=undef){
				if(vdata[j][i][l]>maxV) maxV=vdata[j][i][l];
				if(vdata[j][i][l]<minV) minV=vdata[j][i][l];
				
				sumAve+=vdata[j][i][l]; count++;
			}
			
			if(count!=0){
				sumAve/=count;
				avedata[l]=(float)(sumAve);
			}
			
			maxdata[l]=maxV;
			mindata[l]=minV;
			
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++) if(vdata[j][i][l]!=undef){
				double tmp=(vdata[j][i][l]-sumAve);
				sumVar+=tmp*tmp;
			}
			
			if(count!=0){
				sumVar/=count;
				vardata[l]=(float)(sumVar);
			}
		}
		
		return new Variable[]{ave,var,max,min};
	}
}
