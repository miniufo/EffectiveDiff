//
package pack;

import miniufo.application.basic.DynamicMethodsInCTS;
import miniufo.application.statisticsModel.FilterMethods;
import miniufo.descriptor.DataDescriptor;
import miniufo.diagnosis.CartesianSpatialModel;
import miniufo.diagnosis.DiagnosisFactory;
import miniufo.diagnosis.Range;
import miniufo.diagnosis.Variable;
import miniufo.io.DataIOFactory;
import miniufo.io.DataWrite;

//
public final class Alignment{
	//
	private static final String path="d:/Data/MITgcm/barotropicDG/BetaCartRL/";
	
	//
	public static void main(String[] args){
		cAlignment();
	}
	
	static void cAlignment(){
		DiagnosisFactory df=DiagnosisFactory.parseFile(path+"Stat.cts");
		DataDescriptor dd=df.getDataDescriptor();
		
		Range r=new Range("t(1,500)",dd);
		
		Variable[] vs=df.getVariables(r,"tr1","vor");
		
		for(Variable v:vs) v.replaceUndefData(-9999f);
		
		Variable eta=vs[0];
		Variable  sf=vs[1]; for(int i=0;i<4;i++) FilterMethods.smooth9(sf);
		
		CartesianSpatialModel csm=new CartesianSpatialModel(dd);
		DynamicMethodsInCTS dm=new DynamicMethodsInCTS(csm);
		
		Variable[] grdE=dm.c2DGradient(eta);
		Variable[] grdS=dm.c2DUnitGradient(sf);
		
		Variable align=projection(grdE,grdS); align.setName("align");
		Variable grdmg=grdE[0].hypotenuseEq(grdE[1]); grdmg.setName("grd");
		
		DataWrite dw=DataIOFactory.getDataWrite(dd,path+"Alignment/Align4.dat");
		dw.writeData(dd,eta,sf,align,grdmg); dw.closeFile();
	}
	
	static Variable projection(Variable[] grd1,Variable[] grd2){
		int t=grd1[0].getTCount(); int z=grd1[0].getZCount();
		int y=grd1[0].getYCount(); int x=grd1[0].getXCount();
		
		float undef=grd1[0].getUndef();
		
		Variable align=new Variable("align",grd1[0]);
		align.setCommentAndUnit("alignment (1)");
		align.setValue(undef);
		
		float[][][][] x1data=grd1[0].getData();
		float[][][][] y1data=grd1[1].getData();
		float[][][][] x2data=grd2[0].getData();
		float[][][][] y2data=grd2[1].getData();
		float[][][][] aldata=  align.getData();
		
		if(grd1[0].isTFirst()){
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(x1data[l][k][j][i]!=undef&&y1data[l][k][j][i]!=undef&&x2data[l][k][j][i]!=undef&&y2data[l][k][j][i]!=undef)
				aldata[l][k][j][i]=x1data[l][k][j][i]*y2data[l][k][j][i]-y1data[l][k][j][i]*x2data[l][k][j][i];
			
		}else{
			for(int l=0;l<t;l++)
			for(int k=0;k<z;k++)
			for(int j=0;j<y;j++)
			for(int i=0;i<x;i++)
			if(x1data[k][j][i][l]!=undef&&y1data[k][j][i][l]!=undef&&x2data[k][j][i][l]!=undef&&y2data[k][j][i][l]!=undef)
				aldata[k][j][i][l]=x1data[k][j][i][l]*y2data[k][j][i][l]-y1data[k][j][i][l]*x2data[k][j][i][l];
		}
		
		return align.absEq();
	}
}
