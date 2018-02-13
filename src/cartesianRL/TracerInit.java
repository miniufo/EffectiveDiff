//
package cartesianRL;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteOrder;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import common.MITgcmUtil;
import common.MITgcmUtil.DataPrec;
import miniufo.basic.ArrayUtil;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.CtlDataWriteStream;
import miniufo.io.IOUtil;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.PI;

//
public final class TracerInit{
	//
	private static final int y=Parameters.y;
	private static final int x=Parameters.x;
	
	private static final float depth=Parameters.depth;
	private static final float tauMax=Parameters.tauMax;
	
	private static final String path=Parameters.path;
	
	
	//
	public static void main(String[] args){
		//Variable v=DiagnosisFactory.getVariables(path+"Leith4/DYE/dye_SF.ctl","","dye")[0];
		//replaceTracerFields(path+"Leith4/DYE/pickup_ptracers.0000103680.data",v,DataPrec.float64);
		//generateBath(path+"BATH/bath.dat");
		generateTaux(path+"EXF/taux.dat");
		//for(int i=1;i<=16;i++) generateDye(path+"DYE/dye"+i+".dat",i);
		//generateDyeFromSF(path+"Leith10/meanSF.cts",path+"Leith10/DYE/dye_SF.dat");
	}
	
	static void generateDyeFromSF(String meanSF,String instSF,String out){
		Variable sfm=DiagnosisFactory.getVariables(meanSF,"","sf")[0]; // time mean SF
		Variable sfi=DiagnosisFactory.getVariables(instSF,"","sf")[0]; // instanteous SF
		
		float[][] mdata=sfm.getData()[0][0];
		float[][] idata=sfi.getData()[0][0];
		
		float[] exm=ArrayUtil.getExtrema(mdata);
		float[] exi=ArrayUtil.getExtrema(idata);
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++){
			mdata[j][i]=-1f+2f*(mdata[j][i]-exm[0])/(exm[1]-exm[0]);
			idata[j][i]=-1f+2f*(idata[j][i]-exi[0])/(exi[1]-exi[0]);
		}
		
		sfm.plusEq(2);
		sfi.plusEq(2);
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(out,ByteOrder.BIG_ENDIAN);
		cdws.writeData(sfm,sfi); cdws.closeFile();
	}
	
	static void replaceTracerFields(String pickup,Variable v,DataPrec prec){
		int fields=IOUtil.getFileLength(pickup),fieldLen=x*y*(prec==DataPrec.float32?4:8);
		
		if(fields%fieldLen!=0)
		throw new IllegalArgumentException("invalid file length ("+fields+"), should be multiples of "+fieldLen);
		
		fields/=fieldLen;
		
		try(RandomAccessFile rafR=new RandomAccessFile(pickup,"r" );
			RandomAccessFile rafW=new RandomAccessFile(IOUtil.getCompleteFileNameWithoutExtension(pickup)+".modified","rw")){
			List<float[][]> data=Stream.generate(()->MITgcmUtil.readFloatBE(rafR.getChannel(),x,y,prec)).limit(fields).collect(Collectors.toList());
			
			for(int i=0,I=data.size()-2;i<I;i++) data.set(i,v.getData()[0][0]);
			
			data.stream().forEach(array->MITgcmUtil.writeFloatBE(rafW.getChannel(),array,prec));
			
		}catch(IOException e){ e.printStackTrace(); System.exit(0);}
	}
	
	static void generateBath(String fname){
		Variable bath=new Variable("bath",new Range(1,1,y,x));
		
		float[][] bdata=bath.getData()[0][0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		if(i!=x-1&&j!=y-1) bdata[j][i]=-depth; // set walls
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(fname,ByteOrder.BIG_ENDIAN);
		cdws.writeData(bath); cdws.closeFile();
	}
	
	static void generateDye(String fname,int num){
		Variable dye=new Variable("dye"+num,new Range(1,1,y,x));
		
		float[][] tdata=dye.getData()[0][0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++)
		tdata[j][i]=(float)(1.0000312*sin(2.0*PI*j/(y-2.0))*sin(PI*i/(x-2.0))+1.0+num);
		
		System.out.println("dye "+num+": ["+ArrayUtil.getMin(tdata)+", "+ArrayUtil.getMax(tdata)+"]");
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(fname,ByteOrder.BIG_ENDIAN);
		cdws.writeData(dye); cdws.closeFile();
	}
	
	static void generateTaux(String fname){
		Variable taux=new Variable("taux",new Range(1,1,y,x));
		
		float[][] tdata=taux.getData()[0][0];
		
		for(int j=0;j<y;j++)
		for(int i=0;i<x;i++) tdata[j][i]=(float)((1.00001-cos(2.0*PI*j/(y-2.0)))*tauMax/2.0);
		
		CtlDataWriteStream cdws=new CtlDataWriteStream(fname,ByteOrder.BIG_ENDIAN);
		cdws.writeData(taux); cdws.closeFile();
	}
}
