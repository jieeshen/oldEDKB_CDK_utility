/**
 * Created with IntelliJ IDEA.
 * User: JShen
 * Date: 2/27/13
 * Time: 11:33 PM
 * This program is calculate the CDK_FINGERPRINT and other substructure descriptors used in EDKB and LTKB.
 * It reads a table which contains a column of smiles and write out the following descriptors:
 *      ID	SMILES	CDK_FINGERPRINT	NATOM	NBOND	ATOMLIST	BOND1_ATOM	BOND2_ATOM	BOND_TYPE	TOPOINDEX
 *
 *
 */

import org.openscience.cdk.Molecule;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.io.MDLWriter;
import org.openscience.cdk.smiles.SmilesParser;

import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.fingerprint.Fingerprinter;
import javax.vecmath.*;


import java.io.*;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.StringTokenizer;

//import java.io.PipedReader;

public class oldCdkFP_gen_v2 {
    private static final String[] element = {
            "C",
            "O",
            "H",
            "S",
            "N",
            "Cl",
            "Br",
            "P",
            "I",
            "Cu",
            "Si",
            "F",
            "Se",
            "Hg",
            "B",
            "Sb",
            "Bi",
            "Ti",
            "As",
            "Zn",
            "Ni",
            "Sn",
            "Cr",
            "Cd",
            "Co",
            "Mn",
            "Fe",
            "Al",
            "Na",
            "Ce",
            "Be",
            "Mg",
            "Th",
            "V",
            "Zr",
            "Pt",
            "Te",
            "Pb",
            "Ge",
            "K",
            "Li",
            "Tl",
            "Cs",
            "Pd",
            "Ag",
            "Ca",
            "Ba",
            "U",
            "Mo",
            "Ru",
            "Au",
            "Rb",
            "Pr",
            "Sm",
            "Ga",
            "Rh",
            "Hf",
            "W",
            "Os",
            "Nb",
            "Ta",
            "In",
            "Lu",
            "Sc",
            "Sr",
            "Yb",
            "Tb",
            "Eu",
            "Gd",
            "Dy",
            "Nd",
            "Tm",
            "La",
            "Y",
            "Re",
            "Ir",
            "Er",
            "Ho",
            "Ac",
            "D",
            "Fr",
            "He",
            "Po",
            "Am"
    };

    public static ArrayList<String []> readList(String fileName, String separator) throws Exception {
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fileName)));
        String line = null;
        ArrayList<String[]> outList = new ArrayList<String[]>();
        int i=0;
        while((line=br.readLine())!=null){
            String[] data=line.split(separator);
            outList.add(data);
            i++;
        }
        return outList;
    }

    static	public  String getMolStructure(String m)
    {


        int topoindex=0;

        StringTokenizer mst,st;
        mst=new StringTokenizer (m,"\n");
        mst.nextToken ();
        String s;

        s=mst.nextToken();
//        System.out.println(s);
        int natom,nbond;
        String satom,sbond,atomlist,bond1,bond2,bondtype,stmp;
        //	 st=new StringTokenizer (s," ");
        satom=s.substring (0,3).trim ();
//        System.out.println(satom);
        sbond =s.substring (3,6).trim ();
//        System.out.println(sbond);
        atomlist="";
        bond1="";
        bond2="";
        bondtype="";
        natom=Integer.parseInt (satom);
        nbond=Integer.parseInt (sbond);
        int i; int tmp;
        //System.out.println("error3");
        for(i=0;i<natom;i++)
        {
            s=mst.nextToken  ();

            //	  st=new StringTokenizer (s," ");
            //	  st.nextToken ();
            //	  st.nextToken ();
            //	  st.nextToken ();
            stmp=s.substring (31,34).trim ();
            //System.out.println("read="+stmp);
            tmp=checkAtom(stmp);
            if (tmp <0) {
                System.out.println("Illegal atom!");
                return "error";
                //	  return null;
            }

            if (i==0)
                atomlist=""+tmp;
            else
                atomlist=atomlist+","+tmp;
        }
        for(i=0;i<nbond;i++)
        {
            s=mst.nextToken  ();
            //	   st=new StringTokenizer (s," ");
            if (i==0)
            {
                bond1=s.substring (0,3).trim();

                bond2=s.substring (3,6).trim();
                bondtype=s.substring (6,9).trim();
            }
            else

            {
                bond1=bond1+","+s.substring (0,3).trim();
                bond2=bond2+","+s.substring (3,6).trim();
                bondtype=bondtype+","+s.substring (6,9).trim();
            }
        }

        st=new StringTokenizer (atomlist,",");
        while(st.hasMoreTokens ())
        {
            tmp=Integer.parseInt (st.nextToken ());
            topoindex=topoindex+tmp*tmp;
        }
        st=new StringTokenizer (bondtype,",");
        while(st.hasMoreTokens ())
        {
            tmp=Integer.parseInt (st.nextToken ());
            topoindex=topoindex+tmp*tmp;
        }

		String outputStr=(natom+"\t"+nbond+"\t"+atomlist+"\t"+bond1+"\t"+bond2+"\t"+bondtype+"\t"+topoindex);
//        System.out.print(natom+"\t"+nbond+"\t"+atomlist+"\t"+bond1+"\t"+bond2+"\t"+bondtype+"\t"+topoindex);
        return outputStr;
    }
    static	int checkAtom(String s)
    {
        int i,find=-1;
        boolean exist=false;
        for(i=0;i<element.length;i++)
        { if(s.equals (element[i]))
        {find=i;break;}
        }
        return find;
    }


    public static void main(String[] args) throws Exception{


        String inputFile="C:\\JShen\\Research\\EDKB\\ER\\DEAC\\130227toGerry\\COMPOUND.txt";
        String outputFile="C:\\JShen\\Research\\EDKB\\ER\\DEAC\\130227toGerry\\DEAC_8212_CDK_FP_old.txt";
        ArrayList<String[]> dataMatrix=readList(inputFile,"\t");
        String[] titles=dataMatrix.get(0);
        int smiLoc=0;
        for(int i=0;i<titles.length;i++){
            if (titles[i].equals("SMILES")){
                smiLoc=i;
                break;
            }
        }
        int noInp=dataMatrix.size()-1;

        PrintWriter pw=new PrintWriter(new OutputStreamWriter(new FileOutputStream(outputFile)),true);
        String smi=new String();

        Molecule mm=null;
        SmilesParser sp = new SmilesParser();

        for(int i=1;i<=noInp;i++){
            String[] dataString=dataMatrix.get(i);
            smi=dataString[smiLoc];

//            for (int j=0;j<titles.length;j++){
//                try{
//                    pw.printf("%s\t", dataMatrix.get(i)[j]);
//                }catch (Exception e){
//                    pw.printf("%s\t", "");
//                }
//            }
            pw.printf("%d\t",i);
            pw.printf("%s\t",smi);
            if(smi.length()>250){
                System.out.printf("%d\tTOO LARGE\n",i);
                pw.println("Molecule is too Large");
                continue;
            }
            //System.out.println(i);

            try{
                mm = sp.parseSmiles(smi);
            }catch(Exception e){
                System.out.printf("%d\tCan not handle\n",i);
                pw.println("Can not handle");
                continue;
            }
    //        int nbond=mm.getBondCount();
    //        int natom=mm.getAtomCount();
    //        System.out.println(nbond+" "+natom);

            Fingerprinter fingerprinter = new Fingerprinter();
            try{
                BitSet fingerprint = fingerprinter.getFingerprint(mm);
                int l=fingerprint.length();
                StringBuffer sb=new StringBuffer (l);

                for (int j=0;j<l;j++){
                    if(fingerprint.get(j)){
                        sb.append('1');
                    }else{
                        sb.append('0');
                    }

                }
                pw.print(sb);
                pw.print("\t");
            }catch(Exception e){
                System.out.printf("%d\tCan not generate FP\n",i);
                pw.println("Can not generate FP");
                continue;
            }

            Bond[] bonds = (  Bond[]) mm.getBonds();
            for (int j = 0; j < bonds.length; j++) {
                if (bonds[j].getFlag(CDKConstants.ISAROMATIC))
                    bonds[j].setOrder(4.0);}


            StringWriter sw=new StringWriter ();
            MDLWriter mdlWriter =new MDLWriter(sw);
            mdlWriter.writeMolecule(mm);
            mdlWriter.close();
            String mol=sw.toString ();

            String subDesStr=getMolStructure(mol);
            pw.print(subDesStr);
            pw.println();
        }
    }
}
