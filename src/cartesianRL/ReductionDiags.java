//
package cartesianRL;

import java.util.stream.Stream;
import miniufo.application.basic.DynamicMethodsInCTS;
import miniufo.application.basic.VelocityFieldInCTS;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.CartesianSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;


//
public final class ReductionDiags{
	//
	private static final DiagnosisFactory df=DiagnosisFactory.parseFile("F:/cartRL_advSchemes/Leith1_k0/Stat.cts");
	private static final DataDescriptor dd=df.getDataDescriptor();
	private static final CartesianSpatialModel csm=new CartesianSpatialModel(dd);
	private static final DynamicMethodsInCTS dm=new DynamicMethodsInCTS(csm);
	
	private static final int y=dd.getYCount();
	private static final int x=dd.getXCount();
	
	private static final String path=Grids.path;
	
	
	public static void main(String[] args){
		df.setPrinting(false);
		
		//cMeanAndLastTimeStreamFunction("Leith2/meanSF.dat",301,3600);
		
		/*** get tracer mean, variance, and extrema **
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"Leith1/Leith1_k100/aveVarEx.dat"); cdws.setPrinting(false);
		df.getVariablesTimeByTime("tr1","tr2","tr3","tr4","tr5","tr6","tr7","tr8","tr9","tr10")
			.peek(vs->{
				int t=vs[0].getRange().getTRange()[0];
				if(t%100==0) System.out.println(t);
			})
			.flatMap(vs->reduceToMeanVarianceAndExtremes(vs))
			.forEach(v->cdws.writeData(v));
		cdws.closeFile();*/
		
		/*** get squared tracer and tracer gradient **
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"Leith1/Leith1_k100/squaredGrd.dat"); cdws.setPrinting(false);
		df.getVariablesTimeByTime("tr1","tr2","tr3","tr4","tr5","tr6","tr7","tr8","tr9","tr10")
			.peek(vs->{
				int t=vs[0].getRange().getTRange()[0];
				if(t%100==0) System.out.println(t);
			})
			.flatMap(vs->reduceToSquareAndGrdSquareSum(vs))
			.forEach(v->cdws.writeData(v));
		cdws.closeFile();*/
		
		/*** get kinetic energy ***/
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"Leith1/Leith1_k100/KE.dat"); cdws.setPrinting(false);
		df.getVariablesTimeByTime("u","v")
			.peek(vs->{
				int t=vs[0].getRange().getTRange()[0];
				if(t%100==0) System.out.println(t);
			})
			.map(vs->reduceToKESum(vs[0],vs[1]))
			.forEach(v->cdws.writeData(v));
		cdws.closeFile();
	}
	
	
	static void cMeanAndLastTimeStreamFunction(String out,int tstr,int tend){
		/*** get mean u,v ***/
		Variable um=df.getVariableTimeByTime(tstr,tend,"u").reduce((v1,v2)->v1.plusEq(v2)).get();
		Variable vm=df.getVariableTimeByTime(tstr,tend,"v").reduce((v1,v2)->v1.plusEq(v2)).get();
		
		um.divideEq(tend-tstr+1);
		vm.divideEq(tend-tstr+1);
		
		/*** get last time u,v ***/
		Variable ui=df.getVariables(new Range("t("+tend+","+tend+")",dd),"u")[0];
		Variable vi=df.getVariables(new Range("t("+tend+","+tend+")",dd),"v")[0];
		
		CartesianSpatialModel csm=new CartesianSpatialModel(dd);
		VelocityFieldInCTS     vf=new VelocityFieldInCTS(csm);
		
		Variable sfm=vf.cStreamFunctionBySOR(um,vm);
		Variable sfi=vf.cStreamFunctionBySOR(ui,vi);
		
		um.setName("um"); ui.setName("ui"); sfm.setName("sfm");
		vm.setName("vm"); vi.setName("vi"); sfi.setName("sfi");
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+out);
		dw.writeData(dd,um,vm,sfm,ui,vi,sfi); dw.closeFile();
	}
	
	
	
	static Stream<Variable> reduceToSquareAndGrdSquareSum(Variable[] vs){
		return Stream.of(vs).flatMap(v->reduceToSquareAndGrdSquareSum(v));
	}
	
	static Stream<Variable> reduceToMeanVarianceAndExtremes(Variable[] vs){
		return Stream.of(vs).flatMap(v->reduceToMeanVarianceAndExtremes(v));
	}
	
	static Variable reduceToKESum(Variable u,Variable v){
		changeBCToUndef(u);
		changeBCToUndef(v);
		
		Variable ke=u.square().plusEq(v.square()).divideEq(2f);
		ke.setName("ke"); ke.setComment("kinetic energy");
		
		float undef=ke.getUndef();
		
		float[][] vdata=ke.getData()[0][0];
		
		double sum=0; int count=0;
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) if(vdata[j][i]!=undef){
			sum+=vdata[j][i]; count++;
		}
		
		Variable keRe=new Variable(ke.getName(),new Range(1,1,1,1));
		keRe.getRange().setTRange(v.getRange());
		
		keRe.getData()[0][0][0][0]=(float)sum;
		
		if(count==0) throw new IllegalArgumentException("count is zero");
		
		return keRe;
	}
	
	static Stream<Variable> reduceToSquareAndGrdSquareSum(Variable v){
		changeBCToUndef(v);
		
		Variable vgrd=dm.c2DGradientMagnitude(v);
		
		float undef=v.getUndef();
		
		Variable vari2Sum=new Variable(v.getName()+ "2sum",new Range(1,1,1,1));
		Variable vgrd2Sum=new Variable(v.getName()+"g2sum",new Range(1,1,1,1));
		
		vari2Sum.setCommentAndUnit("sum of squared "+v.getName()+" (unit)"  );
		vgrd2Sum.setCommentAndUnit("sum of squared "+v.getName()+" gradient (unit^2)");
		
		vari2Sum.getRange().setTRange(v.getRange());
		vgrd2Sum.getRange().setTRange(v.getRange());
		
		float[][] vdata=   v.getData()[0][0];
		float[][] gdata=vgrd.getData()[0][0];
		float[] vrdata=vari2Sum.getData()[0][0][0];
		float[] vgdata=vgrd2Sum.getData()[0][0][0];
		
		double sum=0; int count=0;
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) if(vdata[j][i]!=undef){
			sum+=vdata[j][i]*vdata[j][i]; count++;
		}
		
		if(count!=0){
			if(count!=(x-1)*(y-1)) System.out.println("tr count: "+count+" "+v.getName());
			vrdata[0]=(float)sum;
			
		}else throw new IllegalArgumentException("count is 0");
		
		sum=0; count=0;
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) if(gdata[j][i]!=undef){
			sum+=gdata[j][i]*gdata[j][i]; count++;
		}
		
		if(count!=0){
			if(count!=(x-2)*(y-2)) System.out.println("grd count: "+count+" "+v.getName());
			vgdata[0]=(float)sum;
			
		}else throw new IllegalArgumentException("count is 0");
		
		return Stream.of(vari2Sum,vgrd2Sum);
	}
	
	static Stream<Variable> reduceToMeanVarianceAndExtremes(Variable v){
		float undef=v.getUndef();
		
		Variable ave=new Variable(v.getName()+"ave",new Range(1,1,1,1));
		Variable var=new Variable(v.getName()+"var",new Range(1,1,1,1));
		Variable max=new Variable(v.getName()+"max",new Range(1,1,1,1));
		Variable min=new Variable(v.getName()+"min",new Range(1,1,1,1));
		
		ave.setCommentAndUnit("average  of "+v.getName()+" (unit)"  );
		var.setCommentAndUnit("variance of "+v.getName()+" (unit^2)");
		max.setCommentAndUnit("maximum  of "+v.getName()+" (unit)"  );
		min.setCommentAndUnit("minimum  of "+v.getName()+" (unit)"  );
		
		ave.getRange().setTRange(v.getRange());
		var.getRange().setTRange(v.getRange());
		max.getRange().setTRange(v.getRange());
		min.getRange().setTRange(v.getRange());
		
		float[][] vdata=  v.getData()[0][0];
		float[] avedata=ave.getData()[0][0][0];
		float[] vardata=var.getData()[0][0][0];
		float[] maxdata=max.getData()[0][0][0];
		float[] mindata=min.getData()[0][0][0];
		
		float maxV=Float.MIN_VALUE;
		float minV=Float.MAX_VALUE;
		
		double sumVar=0,sumAve=0; int count=0;
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) if(vdata[j][i]!=undef){
			if(vdata[j][i]>maxV) maxV=vdata[j][i];
			if(vdata[j][i]<minV) minV=vdata[j][i];
			
			sumAve+=vdata[j][i]; count++;
		}
		
		if(count!=0){
			if(count!=(x-1)*(y-1)) System.out.println(count+v.getName());
			sumAve/=count;
			avedata[0]=(float)(sumAve);
			
		}else throw new IllegalArgumentException("count is 0");
		
		maxdata[0]=maxV;
		mindata[0]=minV;
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) if(vdata[j][i]!=undef){
			double tmp=(vdata[j][i]-sumAve);
			sumVar+=tmp*tmp;
		}
		
		if(count!=0){
			sumVar/=count;
			vardata[0]=(float)(sumVar);
		}
		
		if(avedata[0]>maxdata[0]) throw new IllegalArgumentException("ave > max");
		if(avedata[0]<mindata[0]) throw new IllegalArgumentException("ave < min");
		
		return Stream.of(ave,var,max,min);
	}
	
	
	static void changeBCToUndef(Variable v){
		int t=v.getTCount(); int z=v.getZCount();
		int y=v.getYCount(); int x=v.getXCount();
		
		float[][][][] vdata=v.getData();
		
		if(v.isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++) vdata[l][k][j][x-1]=-9.99e8f;
				for(int i=0;i<x;i++) vdata[l][k][y-1][i]=-9.99e8f;
			}
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++){
				for(int j=0;j<y;j++) vdata[k][j][x-1][l]=-9.99e8f;
				for(int i=0;i<x;i++) vdata[k][y-1][i][l]=-9.99e8f;
			}
		}
		
		v.setUndef(-9.99e8f);
	}
}
