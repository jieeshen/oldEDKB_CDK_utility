/**
 * Created with IntelliJ IDEA.
 * User: JShen
 * Date: 2/27/13
 * Time: 11:33 PM
 * To change this template use File | Settings | File Templates.
 */

import org.openscience.cdk.Molecule;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.io.*;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor;
import org.openscience.cdk.smiles.*;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.aromaticity.HueckelAromaticityDetector;
import org.openscience.cdk.renderer.Renderer2DModel;
import org.openscience.cdk.io.MDLWriter;
import org.openscience.cdk.smiles.*;
import org.openscience.cdk.exception.*;

import java.awt.event.ActionEvent;
import org.openscience.cdk.Molecule;
import java.awt.event.*;
//import java.io.PipedReader;

import java.util.*;

import org.openscience.cdk.io.MDLWriter;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.config.IsotopeFactory;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.Bond;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.fingerprint.Fingerprinter;


import java.io.*;

public class oldCdkFP_gen {
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

    static	public   void getMolStructure(String m)
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
                return;
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

//		System.out.print(natom+"\t"+nbond+"\t"+atomlist+"\t"+bond1+"\t"+bond2+"\t"+bondtype+"\t"+topoindex+"\n");
        System.out.print(natom+"\t"+nbond+"\t"+atomlist+"\t"+bond1+"\t"+bond2+"\t"+bondtype+"\t"+topoindex);

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
        Molecule mm=null;
        SmilesParser sp = new SmilesParser();

        mm = sp.parseSmiles("CCCCCCCCCc1ccc(O)cc1");
//        int nbond=mm.getBondCount();
//        int natom=mm.getAtomCount();
//        System.out.println(nbond+" "+natom);

        Fingerprinter fingerprinter = new Fingerprinter();
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

        System.out.print(sb);
        System.out.print("\t");

        Bond[] bonds = (  Bond[]) mm.getBonds();
        for (int i = 0; i < bonds.length; i++) {
            if (bonds[i].getFlag(CDKConstants.ISAROMATIC))
                bonds[i].setOrder(4.0);}

        StringWriter sw=new StringWriter ();
        MDLWriter mdlWriter =new MDLWriter(sw);
        mdlWriter.writeMolecule(mm);
        mdlWriter.close();
        String mol=sw.toString ();
        System.out.println(mol);

        getMolStructure(mol);
        System.out.println();

    }
}
