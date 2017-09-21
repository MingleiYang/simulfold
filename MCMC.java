import java.util.*;
import java.io.*;

public class MCMC{

    public Alignment[] alignment;
    public Tree[] tree; // initialized in SimulFold.main()
    public double[] temperature;

    public MCMC(int numberofParalelChains, double maxT){
	alignment = new Alignment[numberofParalelChains];
	tree = new Tree[numberofParalelChains];
	temperature = new double[numberofParalelChains];
	for(int i = 0 ; i < numberofParalelChains; i++){
	    alignment[i] = new Alignment();
	    temperature[i] = Math.pow(Math.pow(maxT,1.0/(numberofParalelChains - 1)),i);
	    System.out.println("temperature: "+temperature[i]);
	}
    }

    public int sampleTreeEdge(int i) {

	alignment[i].deleteValidation();
	double oldLogLikelihood = alignment[i].logLikelihood(tree[i],temperature[i]);
	TreeNode sampledNode = tree[i].selectNode();
	double oldLength = sampledNode.changeEdgeLength(); // sample edge length and return old edge length
	alignment[i].deleteValidation();
	double newLogLikelihood = alignment[i].logLikelihood(tree[i],temperature[i]);
	Random random = new Random();

	if(random.nextDouble() < Math.exp(newLogLikelihood - oldLogLikelihood + sampledNode.edgeLength - oldLength) * 
	   (sampledNode.edgeLength < SimulFold.spanForEdgeSampling ? 
	    sampledNode.edgeLength + SimulFold.spanForEdgeSampling :
	    2.0 * SimulFold.spanForEdgeSampling) /
	   (oldLength < SimulFold.spanForEdgeSampling ?
	    oldLength + SimulFold.spanForEdgeSampling : 
	    2.0 * SimulFold.spanForEdgeSampling)){
	    // acceptance, do nothing
	    System.out.println("Edge Sampling on temperature["+i+"]: ACCEPTED\t"/*+
			       oldLogLikelihood+" "+newLogLikelihood+" "+
			       oldLength+" "+sampledNode.edgeLength*/);
	    return 1;
	    
	}
	else{
	    // rejection, restore the old length
	    alignment[i].deleteValidation();
	    sampledNode.restoreEdgeLength(oldLength);
	    System.out.println("Edge Sampling on temperature["+i+"]: REJECTED\t"/*+
			       oldLogLikelihood+"\t"+newLogLikelihood+" "+
			       oldLength+" "+sampledNode.edgeLength*/);
	    return -1;
	}
    }


    public int nni(int i){ // Nearest Neighbour Interchange

	alignment[i].deleteValidation();
	double oldLogLikelihood = alignment[i].logLikelihood(tree[i],temperature[i]);
	TreeNode sampledNode;

	// select a node which is not the root and its parent is not ther root
	while((sampledNode = tree[i].selectNode()).parent == null || 
	      sampledNode.parent.parent == null);
	
	///////////////////////////
	// swap the aunt and sampled nodes 
	// We could move these lines to the else-section below
	TreeNode motherNode = sampledNode.parent;
	TreeNode grandmotherNode = motherNode.parent;
	TreeNode auntNode = (grandmotherNode.leftChild == motherNode ? 
			     grandmotherNode.rightChild : grandmotherNode.leftChild);

	
	auntNode.parent = motherNode; // change tree
	if(motherNode.leftChild == sampledNode){
	    motherNode.leftChild = auntNode;// change tree
	}
	else{
	    motherNode.rightChild = auntNode;// change tree
	}
	
	sampledNode.parent = grandmotherNode; // change tree
	if(grandmotherNode.leftChild == auntNode){
	    grandmotherNode.leftChild = sampledNode; // change tree
	}
	else{
	    grandmotherNode.rightChild = sampledNode; // change tree
	}
	/////////////////////////////
	alignment[i].deleteValidation();
	double newLogLikelihood = alignment[i].logLikelihood(tree[i],temperature[i]);
	Random random = new Random();	
	if(random.nextDouble() < Math.exp(newLogLikelihood - oldLogLikelihood)){
	    //accept the tree, do nothing
	    System.out.println("NNI Sampling on temperature["+i+"]: ACCEPTED\t"/*+
							   oldLogLikelihood+" "+newLogLikelihood*/);
	    return 1;
	   
	}
	else{
	    //restore the tree
	    auntNode.parent = grandmotherNode;
	    if(grandmotherNode.leftChild == sampledNode){
		grandmotherNode.leftChild = auntNode;
	    }
	    else{
		grandmotherNode.rightChild = auntNode;
	    }
	    sampledNode.parent = motherNode;
	    if(motherNode.leftChild == auntNode){
		motherNode.leftChild = sampledNode;
	    }
	    else{
		motherNode.rightChild = sampledNode;
	    }
	    System.out.println("NNI Sampling on temperature["+i+"]: REJECTED\t"/*+
							   oldLogLikelihood+" "+newLogLikelihood*/);
	    return -1;
	}
    }

    public int sampleAlignment(int i){

	Random generator = new Random();

	alignment[i].deleteValidation();
	double oldLogLikelihood = alignment[i].logLikelihood(tree[i],temperature[i]);

	// Find the possible intervals to be resampled
	HashSet<SingleStranded> intervalSet = new HashSet<>(); // HashSet contains SingleStranded objects
	AlignmentColumn actual = alignment[i].first;

	while(actual != null){

	    SingleStranded interval = new SingleStranded();
	    while(actual != null && actual.basepair == null){
		if(interval.first == null){
		    interval.first = actual; // set values of the SingleStranded object
		}
		interval.last = actual; // set values of the SingleStranded object
		interval.length++;      // set values of the SingleStranded object
		actual = actual.next;

		
	    }
	    if(interval.length != 0){
		intervalSet.add(interval);
	    }
	    if(actual != null){
		actual =actual.next;
	    }
	}	
	// choose an interval
	int choosenIndex = generator.nextInt(intervalSet.size());
	// draw an object unifromly form a HashSet
	Iterator iter = intervalSet.iterator();
	SingleStranded choosenInterval = (SingleStranded) iter.next();
	for(int j = 0; j < choosenIndex; j++){
	    choosenInterval = (SingleStranded) iter.next();
	}

	// creating a new alignment
	AlignmentProposal proposal = tree[i].root.alignmentProposal(alignment[i], choosenInterval.first, choosenInterval.last);
	proposal.alignment.sequenceName = alignment[i].sequenceName;
            
	// getting the correct permutation of the columns 
	actual = proposal.alignment.first;
	while(actual != null){
	    int[] temp = actual.copyColumn();
	    for(int j = 0; j < temp.length; j++){
		actual.column[proposal.alignment.indexes[j]] = temp[j];
	    }
	    actual = actual.next;
	}
	for(int j = 0; j < proposal.alignment.indexes.length; j++){
	    proposal.alignment.indexes[j] = j;
	}
	proposal.alignment.leaf = new TreeNode[alignment[i].sequenceName.length];
	proposal.alignment.findLeaves(tree[i]);

	//System.out.println("Going to backproposal  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	AlignmentProposal backproposal =
		    tree[i].root.alignmentBackProposal(proposal.alignment, proposal.alignment.first, proposal.alignment.last);
	double checkProb = backproposal.logProb;
	backproposal = 
	    tree[i].root.alignmentBackProposal(alignment[i], choosenInterval.first, choosenInterval.last);

	// putting the resampled alignment back to the alignment
	if(choosenInterval.first.prev != null){
	    choosenInterval.first.prev.next = proposal.alignment.first;
	    proposal.alignment.first.prev = choosenInterval.first.prev;
	}
	else{
	    alignment[i].first = proposal.alignment.first;
	}
	if(choosenInterval.last.next != null){
	    choosenInterval.last.next.prev = proposal.alignment.last;
	    proposal.alignment.last.next = choosenInterval.last.next;
	}
	else{
	    alignment[i].last = proposal.alignment.last;
	}

	/*
	TreeNode[] tempTreeNodes = alignment.leaf;
	alignment.leaf = null;
	alignment.leaf = proposal.alignment.leaf;

	int[] tempIndexes = alignment.indexes;
	alignment.indexes = null;
	alignment.indexes = proposal.alignment.indexes;
	*/


	alignment[i].deleteValidation();
	double newLogLikelihood = alignment[i].logLikelihood(tree[i],temperature[i]);
	if(generator.nextDouble() < Math.exp(newLogLikelihood + backproposal.logProb -
					     (oldLogLikelihood + proposal.logProb))){
	    // Acceptance, do nothing
	    System.out.println("alignment Sampling on temperature["+i+"]: ACCEPTED\t"/*+
			       oldLogLikelihood+" "+newLogLikelihood+" "+
			       backproposal.logProb+" "+proposal.logProb+" "+checkProb*/);
	    return 1;
	}
	else{
	    // Restore the alignment
	    System.out.println("alignment Sampling on temperature["+i+"]: REJECTED\t"/*+
			       oldLogLikelihood+" "+newLogLikelihood+" "+
			       backproposal.logProb+" "+proposal.logProb+" "+checkProb*/);
	    // putting the old sub-alignment back to the alignment
	    if(choosenInterval.first.prev != null){
		choosenInterval.first.prev.next = choosenInterval.first;
	    }
	    else{
		alignment[i].first = choosenInterval.first;
	    }
	    if(choosenInterval.last.next != null){
		choosenInterval.last.next.prev = choosenInterval.last;
	    }
	    else{
		alignment[i].last = choosenInterval.last;
	    }
	    return -1;

	}
	/*
	proposal.alignment.print();
	System.out.println("----------------------------------------------------");
	alignment.print();
	*/
    }

    
    public int sampleStructure(int i){

	double oldStructurePrior = alignment[i].structurePrior(temperature[i]);
	alignment[i].deleteValidation();
	double oldLogLikelihood = alignment[i].logLikelihood(tree[i],temperature[i]);

	// remove some helices
	double removingMH = alignment[i].removeHelices(tree[i],temperature[i]);
	//System.out.println("Structure after removing helices");
	//alignment.first.printStructurefromHere();

	// propose some helices
	double proposingMH = alignment[i].proposeHelices(tree[i],temperature[i]);
	//System.out.println("Structure after proposing helices");
	//alignment.first.printStructurefromHere();

	double newStructurePrior = alignment[i].structurePrior(temperature[i]);
	alignment[i].deleteValidation();
	double newLogLikelihood = alignment[i].logLikelihood(tree[i],temperature[i]);

	Random generator = new Random(); 
	/*
	System.out.println("Metropolis-Hastings ratio: "+(Math.exp(newLogLikelihood - oldLogLikelihood) *
							  removingMH * proposingMH *
							  newStructurePrior / oldStructurePrior));
	*/
	if(generator.nextDouble() < (Math.exp(newLogLikelihood - oldLogLikelihood) *
				     removingMH * proposingMH *
				     newStructurePrior / oldStructurePrior)){
	    // accept the new structure
	    System.out.println("Structure Sampling on temperature["+i+"]: ACCEPTED\t"/*+
			       removingMH+" "+proposingMH+" "+
			       oldLogLikelihood+" "+newLogLikelihood+" "+
			       oldStructurePrior+" "+newStructurePrior*/);
	    /*
	    String[] printed = alignment.printedalignment();
	    for(int j = 0; j < printed.length; j++){
		System.out.println("Alignment\t"+printed[j]);
	    }
	    System.out.println("Structure\t"+alignment.printedcommonstructure());
	    printed = alignment.printedstructure();
	    for(int j = 0; j < printed.length; j++){
		System.out.println("Structure\t"+printed[j]);
	    }
	    */					       

	    return 1;
	}
	else{
	    // reject the new structure, restore the old structure
	    System.out.println("Structure Sampling on temperature["+i+"]: REJECTED\t"/*+
			       removingMH+" "+proposingMH+" "+
			       oldLogLikelihood+" "+newLogLikelihood+" "+
			       oldStructurePrior+" "+newStructurePrior*/);
	    /*
	    String[] printed = alignment.printedalignment();
	    for(int j = 0; j < printed.length; j++){
		System.out.println("Alignment\t"+printed[j]);
	    }
	    System.out.println("Structure\t"+alignment.printedcommonstructure());
	    printed = alignment.printedstructure();
	    for(int j = 0; j < printed.length; j++){
		System.out.println("Structure\t"+printed[j]);
	    }
	    */
	    AlignmentColumn actual = alignment[i].first;
	    while(actual != null){
		actual.basepair = actual.oldBasepair;
		actual.isBasepairAhead = actual.oldIsBasepairAhead;
		actual = actual.next;
	    }
	    /*
	    System.out.println("Structure after restoring helices");
	    printed = alignment.printedalignment();
	    for(int j = 0; j < printed.length; j++){
		System.out.println("Alignment\t"+printed[j]);
	    }
	    System.out.println("Structure\t"+alignment.printedcommonstructure());
	    printed = alignment.printedstructure();
	    */
	    return -1;
	}
    }


    public int swap(int i){
	
	alignment[i].deleteValidation();
	alignment[i+1].deleteValidation();
	double oldLogLikelihood1 = alignment[i].logLikelihood(tree[i],temperature[i]);
	double oldLogLikelihood2 = alignment[i+1].logLikelihood(tree[i+1],temperature[i+1]);
	double oldStructurePrior1 = alignment[i].structurePrior(temperature[i]);
	double oldStructurePrior2 = alignment[i+1].structurePrior(temperature[i+1]);

	alignment[i].deleteValidation();
	alignment[i+1].deleteValidation();
	double newLogLikelihood1 = alignment[i].logLikelihood(tree[i],temperature[i+1]);
	double newLogLikelihood2 = alignment[i+1].logLikelihood(tree[i+1],temperature[i]);
	double newStructurePrior1 = alignment[i].structurePrior(temperature[i+1]);
	double newStructurePrior2 = alignment[i+1].structurePrior(temperature[i]);

	Random generator = new Random();
	if(generator.nextDouble() < (Math.exp(newLogLikelihood1 + newLogLikelihood2 - oldLogLikelihood1 - oldLogLikelihood2) *
				     newStructurePrior1 * newStructurePrior2 / (oldStructurePrior1 * oldStructurePrior2))){
	    //swap things
	    Alignment tempAlignment = alignment[i];
	    alignment[i] = null;
	    alignment[i] = alignment[i+1];
	    alignment[i+1] = null;
	    alignment[i+1] = tempAlignment;

	    Tree tempTree = tree[i];
	    tree[i] = null;
	    tree[i] = tree[i+1];
	    tree[i+1] = null;
	    tree[i+1] = tempTree;
	    System.out.println("Swapping temperature["+i+"] and temperature["+(i+1)+"]: ACCEPTED\t");

	    return 1;
	}
	else{
	    System.out.println("Swapping temperature["+i+"] and temperature["+(i+1)+"]: REJECTED\t");
	    return -1;
	}

    }

    public int[][] completeSampling(){

	// Irmi: rewrite later on

	int[][] acceptance = new int[temperature.length][5];

	// Metropolis-Hastings on the chains
	for(int i = 0; i < temperature.length; i++){
	    if(SimulFold.sampleTree){
		acceptance[i][0] = sampleTreeEdge(i);
		acceptance[i][1] = nni(i);
	    }
	    if(SimulFold.sampleAlignment){
		acceptance[i][2] = sampleAlignment(i);
	    }
	    if(SimulFold.sampleStructure){
		acceptance[i][3] = sampleStructure(i);
	    }
	}

	//Paralel tempering
	
	for(int i = 0; i < temperature.length - 1; i++){
	    acceptance[i][4] = swap(i);
	}
	
	return acceptance;

    }
    
    public void doMCMC(int burnin, int samples, int iteration) throws IOException{

	FileWriter f = new FileWriter("output.txt",false);

	int[][] acceptance;

	//this is for large alignments
	int x = 0;
	while(x < 100*temperature.length){
	    for(int i = 0; i < temperature.length; i++){
		x += sampleStructure(i) + 1;
	    }
	}

	for(int i = 0; i < burnin; i++){
	    acceptance = completeSampling();
	}

	for(int i = 0; i < samples; i++){
	    int[] edge = new int[temperature.length];
	    int[] topology = new int[temperature.length];
	    int[] ali = new int[temperature.length];
	    int[] structure = new int[temperature.length];
	    int[] swap = new int[temperature.length];

	    int[] edgea = new int[temperature.length];
	    int[] topologya = new int[temperature.length];
	    int[] alia = new int[temperature.length];
	    int[] structurea = new int[temperature.length];
	    int[] swapa = new int[temperature.length];
	    for(int j = 0; j < iteration; j++){
		acceptance = completeSampling();
		for(int t = 0; t < temperature.length; t++){
		    if(acceptance[t][0] != 0){
			edge[t]++;
			if(acceptance[t][0] == 1){
			    edgea[t]++;
			}
		    }
		    if(acceptance[t][1] != 0){
			topology[t]++;
			if(acceptance[t][1] == 1){
			    topologya[t]++;
			}
		    }
		    if(acceptance[t][2] != 0){
			ali[t]++;
			if(acceptance[t][2] == 1){
			    alia[t]++;
			}
		    }
		    if(acceptance[t][3] != 0){
			structure[t]++;
			if(acceptance[t][3] == 1){
			    structurea[t]++;
			}
		    }
		}//t, M-H aceptance
		for(int t = 0; t < temperature.length - 1; t++){
		    if(acceptance[t][4] != 0){
			swap[t]++;
			if(acceptance[t][4] == 1){
			    swapa[t]++;
			}
		    }
		}
	    }
	    
	    //for not deviding by 0
	    swap[temperature.length - 1] = 1;

	    double[][] a = new double[temperature.length][5];
	    for(int t = 0; t < temperature.length; t++){
		a[t][0] = (double)edgea[t]/(double)edge[t];
		a[t][1] = (double)topologya[t]/(double)topology[t];
		a[t][2] = (double)alia[t]/(double)ali[t];
		a[t][3] = (double)structurea[t]/(double)structure[t];
		a[t][4] = (double)swapa[t]/(double)swap[t];
	    }

	    report(f,i,a);
	    //System.out.println(mcmc.tree.printedTree());
	}
	f.close();
    }

    public void report(FileWriter f, int reportnumber, double[][] a) throws IOException{
	for(int i = 0; i < temperature.length; i++){
	    alignment[i].deleteValidation();
	    
	    f.write("Temperature["+i+"]\tDiagnosis\tLoglikelihood\t"+alignment[i].logLikelihood(tree[i],temperature[i])+
		    "\tStructure prior (unnormalized)\t"+alignment[i].structurePrior(temperature[i])+
		    "\tEdgeSampling(%)\t"+(a[i][0]*100.0)+
		    "\tTopologySampling(%)\t"+(a[i][1]*100.0)+
		    "\tAlignmentSampling(%)\t"+(a[i][2]*100.0)+
		    "\tStructureSampling(%)\t"+(a[i][3]*100.0)+"\tSwapping(%)\t"+(a[i][4]*100.0)+"\n");
	    f.write("Temperature["+i+"]\tNewick tree\t"+tree[i].printedTree()+"\n");
	    String[] printed = alignment[i].printedalignment();
	    for(int j = 0; j < printed.length; j++){
		f.write("Temperature["+i+"]\tAlignment\t"+printed[j]+"\n");
	    }
	    f.write("Temperature["+i+"]\tAlignment\t\n");
	    f.write("Temperature["+i+"]\tStructure\t"+alignment[i].printedcommonstructure()+"\n");
	    printed = alignment[i].printedstructure();
	    for(int j = 0; j < printed.length; j++){
		f.write("Temperature["+i+"]\tStructure\t"+printed[j]+"\n");
	    }
	}
    }

    public static void main(String[] args) throws IOException{ 
	// run this class with "java MCMC input.fasta tree.newick"
	// for testing the class.
	
	//out of use....
	/*
	MCMC mcmc = new MCMC(5,7.0);
	SimulFold.sequences = mcmc.alignment[0].readAlignment(args[0]);
	mcmc.tree = new Tree(args[1]);
	mcmc.alignment.findLeaves(mcmc.tree);
	mcmc.doMCMC(100,100,10);
	*/
    }

}
