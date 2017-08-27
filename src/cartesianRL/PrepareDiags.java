//
package cartesianRL;

import java.util.Optional;
import java.util.stream.IntStream;
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
public final class PrepareDiags{
	//
	private static final DiagnosisFactory df=DiagnosisFactory.parseFile("I:/Leith2/Stat.cts");
	private static final DataDescriptor dd=df.getDataDescriptor();
	
	private static final int y=dd.getYCount();
	private static final int x=dd.getXCount();
	
	private static final String path=Grids.path;
	
	
	public static void main(String[] args){
		cMeanAndLastTimeStreamFunction("Leith2/meanSF.dat",301,3600);
		
		/*** get tracer mean, variance, and extrema **
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"Leith4/aveVarEx.dat"); cdws.setPrinting(false);
		IntStream.range(1,dd.getTCount()+1).sequential()
			.peek(i->{ if(i%100==0) System.out.println(i);})
			.mapToObj(i->getDataAt(i,"tr1","tr2","tr3","tr4","tr5","tr6","tr7","tr8","tr9","tr10","tr11","tr12","tr13","tr14","tr15","tr16"))
			.map(vs->reduceToMeanVarianceAndExtremes(vs))
			.forEach(vs->vs.forEach(v->cdws.writeData(v)));
		cdws.closeFile();*/
		
		/*** get squared tracer and tracer gradient ***/
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"Leith4/squaredGrd.dat"); cdws.setPrinting(false);
		IntStream.range(1,dd.getTCount()+1).sequential()
			.peek(i->{ if(i%100==0) System.out.println(i);})
			.mapToObj(i->getDataAt(i,"tr1","tr2","tr3","tr4","tr5","tr6","tr7","tr8","tr9","tr10","tr11","tr12","tr13","tr14","tr15","tr16"))
			.map(vs->reduceToSquareAndGrdSquareSum(vs))
			.forEach(vs->vs.forEach(v->cdws.writeData(v)));
		cdws.closeFile();
		
		/*** get kinetic energy **
		CtlDataWriteStream cdws=new CtlDataWriteStream(path+"Leith4/KE.dat"); cdws.setPrinting(false);
		IntStream.range(1,dd.getTCount()+1).sequential()
			.peek(i->{ if(i%100==0) System.out.println(i);})
			.mapToObj(i->getDataAt(i,"u","v"))
			.map(vs->reduceToKineticEnergy(vs))
			.forEach(vs->vs.forEach(v->cdws.writeData(v)));
		cdws.closeFile();*/
	}
	
	
	static void cMeanAndLastTimeStreamFunction(String out,int tstr,int tend){
		/*** get mean u,v ***/
		Variable um=IntStream.range(tstr,tend+1).mapToObj(i->getOneVar(i,"u")).reduce((v1,v2)->v1.plusEq(v2)).get();
		Variable vm=IntStream.range(tstr,tend+1).mapToObj(i->getOneVar(i,"v")).reduce((v1,v2)->v1.plusEq(v2)).get();
		
		um.divideEq(tend-tstr+1);
		vm.divideEq(tend-tstr+1);
		
		/*** get last time u,v ***/
		Variable ui=getOneVar(tend,"u");
		Variable vi=getOneVar(tend,"v");
		
		CartesianSpatialModel csm=new CartesianSpatialModel(dd);
		VelocityFieldInCTS     vf=new VelocityFieldInCTS(csm);
		
		Variable sfm=vf.cStreamFunctionBySOR(um,vm);
		Variable sfi=vf.cStreamFunctionBySOR(ui,vi);
		
		um.setName("um"); ui.setName("ui"); sfm.setName("sfm");
		vm.setName("vm"); vi.setName("vi"); sfi.setName("sfi");
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+out);
		dw.writeData(dd,um,vm,sfm,ui,vi,sfi); dw.closeFile();
	}
	
	
	static Stream<Variable> reduceToKineticEnergy(Stream<Variable> vs){
		Optional<Variable> op=vs.map(v->v.square().divideEq(2f)).reduce((v1,v2)->v1.plusEq(v2));
		
		if(op.isPresent()){
			Variable ke=op.get();
			
			ke.setName("ke");
			ke.setComment("kinetic energy");
			
			return Stream.of(reduceToSum(ke));
			
		}else throw new IllegalArgumentException("no ke found");
	}
	
	static Stream<Variable> reduceToSquareAndGrdSquareSum(Stream<Variable> vs){
		CartesianSpatialModel csm=new CartesianSpatialModel(dd);
		DynamicMethodsInCTS    dm=new DynamicMethodsInCTS(csm);
		
		return vs.flatMap(v->reduceToSquareAndGrdSquareSum(v,dm.c2DGradientMagnitude(v)));
	}
	
	static Stream<Variable> reduceToMeanVarianceAndExtremes(Stream<Variable> vs){
		return vs.flatMap(v->reduceToMeanVarianceAndExtremes(v));
	}
	
	
	static Variable reduceToSum(Variable ke){
		float undef=ke.getUndef();
		
		float[][] vdata=ke.getData()[0][0];
		
		double sum=0; int count=0;
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) if(vdata[j][i]!=undef){
			sum+=vdata[j][i]; count++;
		}
		
		Variable keRe=new Variable(ke.getName(),new Range(1,1,1,1));
		
		keRe.getData()[0][0][0][0]=(float)sum;
		
		if(count==0) throw new IllegalArgumentException("count is zero");
		
		return keRe;
	}
	
	static Stream<Variable> reduceToSquareAndGrdSquareSum(Variable v,Variable vgrd){
		float undef=v.getUndef();
		
		Variable vari2Sum=new Variable(v.getName()+ "2sum",new Range(1,1,1,1));
		Variable vgrd2Sum=new Variable(v.getName()+"g2sum",new Range(1,1,1,1));
		
		vari2Sum.setCommentAndUnit("sum of squared "+v.getName()+" (unit)"  );
		vgrd2Sum.setCommentAndUnit("sum of squared "+v.getName()+" gradient (unit^2)");
		
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
	
	
	static Stream<Variable> getDataAt(int t,String... vars){
		df.setPrinting(false);
		
		Variable[] vs=df.getVariables(new Range("t("+t+","+t+")",dd),vars);
		
		for(Variable v:vs) v.replaceUndefData(-9.99e8f);
		
		return Stream.of(vs);
	}
	
	static Variable getOneVar(int t,String vname){
		df.setPrinting(false);
		
		Variable v=df.getVariables(new Range("t("+t+","+t+")",dd),vname)[0];
		
		v.replaceUndefData(-9.99e8f);
		
		return v;
	}
}
