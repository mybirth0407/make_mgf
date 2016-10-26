package modi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.TreeSet;


public class Spectrum extends ArrayList<Peak> implements Comparable<Spectrum> {

    private String        name = "";    // id
    private double         precursor;
    private int         charge;
    private double         observedMW;    // not used now
    private double         correctedMW; //이온화되기 전 peptide의 MW
    private double        TIC= 0;
    private double        RTTime= 0;
    private double        bpIntensity= 0;
    private double        randomPeakIntensity= 0;
    //2015.9.18 새롭게 추가 ( fast LB 할때 평균, 표준편차 구하기위해서)
    private double        sumOfIntensity = 0;
    private double        sumOfSquaredIntensity = 0;
    private double		  sumOfLastBinIntensity = 0;
    private double 		  sumOfLastBinSquaredIntensity = 0;
    
    private    ArrayList<Peak>     selectedPeak = null;

    public    Spectrum(double precursor, int charge) 
    {
        correctedMW = observedMW = (precursor - Constants.Proton)*charge;
        this.precursor = precursor;
        this.charge = charge;
    }
    public     Spectrum(double precursor, int charge, String name)
    {
        this(precursor, charge);
        this.name = name;
    }
    public    void     setExtraInformation(double bpi, double t){
        bpIntensity= bpi;
        TIC= t;
    }
    public    void     addPeak(Peak p)         { this.add(p); }
    public    double    getObservedMW()    { return observedMW; }
    public    double    getCorrectedMW()    { return correctedMW; }
    public    void    setCorrectedParentMW(double mw)    { 
        correctedMW = mw; 
    }
    //2015.9.18 새롭게 추가
    public double getSumOfIntensity()        { return sumOfIntensity; }
    public double getSumOfSquaredIntensity(){ return sumOfSquaredIntensity; }
    public double getSumOfLastBinIntensity(){ return sumOfLastBinIntensity; }
    public double getSumOfLastBinSquaredIntensity(){ return sumOfLastBinSquaredIntensity;}
    
    public void setSumOfIntensity(double value)            { this.sumOfIntensity = value; }
    public void setSumOfSquaredIntensity(double value)    { this.sumOfSquaredIntensity = value; }
    public void setSumOfLastBinIntensity(double value)		{this.sumOfLastBinIntensity = value; }
    public void setSumOfLastBinSquaredIntensity(double value) {this.sumOfLastBinSquaredIntensity = value;}

    public    double    getPrecursor()            { return precursor; }
    public    double    getBasePaekIntensity()    { return bpIntensity; }
    public    int        getCharge()             { return charge; }
    public    String    getName()                { return name; }
    public    double    getRandomPeakIntensity() { return randomPeakIntensity; }
    public  double  getTIC(){ return TIC; }
    public  double  getRTTime(){ return RTTime; }
    public    Peak    getPeakByMass(double mass)
    {
        Peak result = new Peak(-1, -1, 0);
        for(Peak p : this)
        {
            if(Constants.fEqual(mass, p.getMass()) && p.getIntensity() > result.getIntensity())
                result = p;
        }

        return result;
    }

    public    ArrayList<Peak>        getSelectedPeaks()    { return selectedPeak; }

    public void normalizeIntensityLocally(){
        for(Peak p : this){           
            p.setNormIntensity(p.getIntensity()/getLocalMaxPeak(p.getMass()));    // NormIntensity 는 Peak 의 y 축 값(intensity)에서, 해당 Window 에서 
                                                                                // intensity 가 가장 큰 Peak 의 intensity(y 축) 을 나눈 값이다. 
        }
    }

    /*--- Window 에서 intensity 가 가장 큰 Peak 을 반환 ---*/
    public double getLocalMaxPeak(double center){    // 전달된 인자 center 는 해당 Peak 의 X 축에 해당하는 Mass 값 
        double lMax =0;
        int bin = 50;    // Window 의 크기 
        int index = binarySearchforPeaks( center-bin );    // binarySearchforPeaks 함수에 mass 값에서 Window 의 크기를 뺀값을 인자로 전달하여 index 값설정 

        while(this.get(index).getMass() <= center+bin ){    // index 를 증가시켜주며, 해당 index 의 mass 값이 center + bin 보다 작거나 같으면 계속 실행
            if( this.get(index).getIntensity() > lMax )    // 해당 index 의 y 축이 lMax 값 보다 크다면
                lMax= this.get(index).getIntensity();    // lMax 값은 해당 index 의 y 축 으로 변경 
            index++;    // index 값 증가 
            if( index == this.size() )    // 현재 index 가 Window 의 마지막 Peak 에 해당할때 
                break;
        }   
        return lMax;    // 현재 Window 에서 intensity 가 가장 큰 Peak 반환 
    }

    public double forNewNormalization_getLocalMaxPeak(double center){    // 전달된 인자 center 는 해당 Peak 의 X 축에 해당하는 Mass 값 
        double lMax =0;
        int bin = 25;    // Window 의 크기 
        int index = binarySearchforPeaks( center-bin );    // binarySearchforPeaks 함수에 mass 값에서 Window 의 크기를 뺀값을 인자로 전달하여 index 값설정 

        while( this.get(index).getMass() > center-bin && this.get(index).getMass() <= center+bin ){    // index 를 증가시켜주며, 해당 index 의 mass 값이 center + bin 보다 작거나 같으면 계속 실행
            if( this.get(index).getIntensity() > lMax )    // 해당 index 의 y 축이 lMax 값 보다 크다면
                lMax= this.get(index).getIntensity();    // lMax 값은 해당 index 의 y 축 으로 변경 
            index++;    // index 값 증가 
            if( index == this.size() )    // 현재 index 가 Window 의 마지막 Peak 에 해당할때 
                break;
        }   
        return lMax;    // 현재 Window 에서 intensity 가 가장 큰 Peak 반환 
    }

    public     String    toString() 
    {
        StringBuffer buffer = new StringBuffer();
        buffer.append("Spectrum("+precursor+","+charge+")\n");
        Iterator e = this.iterator();
        while(e.hasNext())
        {
            buffer.append(((Peak)e.next()).toString());
//            buffer.append("\n");
        } 
        return buffer.toString();
    } //global selected size를 observed MW 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    public int peakSelection (double selectionWindowSize, int numOfPeaksInWindow ){
        int globalSelectedSize = (int) (this.observedMW / selectionWindowSize * numOfPeaksInWindow);
        //globalSelectedSize = (int) (this.observedMW / 100 * 2);
        return peakSelection(globalSelectedSize, selectionWindowSize, numOfPeaksInWindow);
    }
    ////globalSelectedSize는 처음 에 나누면, 나눴을때, 윈도우가 몇개가 있는지를 나타낸다. 곱하기 peak의 갯수를 하면, 전체 peak의 갯수
    ////전체 peak의 갯수가 사이즈보다크면 다때려박고 끝나는 경우다. 그러나 만약에 작으면, 그냥 다쓰는 것이다.
    ////우선 백개하고, 제일 작은거 비교해서, 쭉 끝까지 다돈다. global는 인텐시티 크기를 매겼을때, globalselectedSize갯수만큼 높은 peak을 뽑은 것이다. 
    public int peakSelection (int globalSelectedSize, double selectionWindowSize, int numOfPeaksInWindow ){
        assert(this.checkSorted());    //파라미터가 null 값인지 체크
        if(globalSelectedSize > this.size())    //계산된 numOfPeaksInWindow구간 만큼의 peak 의 갯수가, 서브클래스 list 갯수보다 클경우
            globalSelectedSize = this.size();    //globalSelectedSize에 서버클래스 list 갯수값을 넣는다. 
        if(this.size() == 0)    //파라미터의 값이 0일경우(값이 0이 될수는 없음)
            return -1;    //에러 코드
        //일정 비율을 global과 local에 대해서 어떻게 되느냐. 지우는 peak들을 noise라고 한다.
        //-------------------------------global selection부분 (Intensity는 y축 값이다.)-------------------------------// 
        Iterator<Peak> it = this.iterator();    //서브클래스의 Collection 을 Peak 자료구조의 객체로, Iterator Collection 에 넣는다.(변수 it)
        TreeSet<Peak> globalSelected = new TreeSet<Peak>(new IntensityComparator());    //(Intensity 비교)자료구조 Peak 을 TreeSet 클래스를 이용하여 
                                                                                        //오름차순으로 정렬하여 globalSelected에 넣는다. 
        for(int i=0; it.hasNext() && i<globalSelectedSize; i++)    //서브클래스의 값부터 globalSelectedSize만큼
            globalSelected.add(it.next());    //정렬된 TreeSet에 add 
        TreeSet<Peak> notSelected = new TreeSet<Peak>(new MassComparator());    //(mass 비교)선택되지않은 자료구조 Peak 은 notSelected에 넣는다. 
        Peak curPeak;    //임시 Peak형 변수 curPeak(현재 값을 나타내기 위해서) 
        while(it.hasNext())    //it 에 값이 존재할때까지 (즉 선택 안되고 남은놈들)
        {
            curPeak = it.next();    //그다음 peak 을 curPeak에 임시로 넣는다.
            if(curPeak.getIntensity() > globalSelected.first().getIntensity()) // 현재 peak 에서의 intensity 값이, 가장 작은 Intensity 값을 가진 
                                                                               //Peak 보다 클경우
            {
                notSelected.add(globalSelected.first());    //globalSelected의 첫번째 Peak 값을 notSelected에 더해준다. 
                globalSelected.remove(globalSelected.first());    // 가장 작은 peak 을 제거한다. 
                globalSelected.add(curPeak);    //현재 가리키고있는 Peak 을 globalSelected Collection 에 더해준다. 
            }
            else
                notSelected.add(curPeak);    //현재 가리키는 Peak 을 notSelected에 더해준다. 
        }
        //-------------------------------local selection 부분-------------------------------// 
        TreeSet<Peak> selected = new TreeSet<Peak> (new MassComparator());    //(mass 값 비교)자료구조 Peak 을 정렬하여 selected 에 넣는다. 
        selected.addAll(globalSelected);    //global selection 된 모든 List 를 selected 에 추가 
        int currentBinPeakSize;    //현재 Window 의 크기를 나타내기위하여, 변수 설정 
        double currentBinStart, currentBinEnd;    //위치 변수로 사용하기위하여, 변수 설정
        // 추가적인 peak 은 selection 하지 않는다.  in (0, selectionWindowSize)
        for(int i=2; i<(int)(correctedMW/selectionWindowSize+1)*2; i++)   
        {    //한Window를 Bin 이라고 한다. 
            currentBinStart = i*selectionWindowSize/2;    //local selection 할 시작 부분 
            currentBinEnd = currentBinStart + selectionWindowSize;    //local selection 할 끝 부분 
            currentBinPeakSize = selected.subSet(new Peak(0, currentBinStart, 0), new Peak(0, currentBinEnd, 0)).size();
            //Peak 에서 local selection 할 부분의 시작 부분에서, 끝 부분까지 부분의 Peak 의 개수를 subSet 메소드를 이용해서 가져온다. 
            if( currentBinPeakSize < numOfPeaksInWindow ){    //현재 Window 안에 있는 peak 의 갯수 < 한윈도우에 들어가는 peak 의 개수
                                                            //예) 현재 3개뽑아야 하는데, 1개밖에없을 경우 충족시켜줘야한다.
                Peak [] newPeaks = notSelected.subSet(new Peak(0, currentBinStart, 0.), new Peak(0, currentBinEnd, 0.)).toArray(new Peak[0]);
                //notSelected에서 해당 범위의 Peak 을 새로운 newPeaks에 넣는다. 
                Arrays.sort(newPeaks, Collections.reverseOrder(new IntensityComparator()));    //Intensity 순으로, 역순으로 정렬 
                for(int j=0; j<newPeaks.length && j < numOfPeaksInWindow - currentBinPeakSize; j++)   
                {                                    //Bin
                    if( newPeaks[j].getNormIntensity() > Constants.minNormIntensity ){
                        notSelected.remove(newPeaks[j]);    //notSelected에서 newPeaks의 j 에 해당하는 peak 을 제거 
                        selected.add(newPeaks[j]);    //newPeaks의 j 에 해당하는 값을 selected 에 추가 
                    }
                    else break;
                } 
            }
        }
        //Property : b이온/y이온 이냐 
        selected.add(new Peak(-1, Constants.B_ION_OFFSET+Constants.NTERM_FIX_MOD, 0, 1, PeakProperty.N_TERM_B_ION_ONLY));
        selected.add(new Peak(-1, Constants.Y_ION_OFFSET+Constants.CTERM_FIX_MOD, 0, 1, PeakProperty.C_TERM_Y_ION_ONLY));
        selected.add(new Peak(-1, correctedMW-Constants.H2O+Constants.Proton-Constants.CTERM_FIX_MOD, 0, 1, PeakProperty.C_TERM_B_ION_ONLY));
        selected.add(new Peak(-1, correctedMW+Constants.Proton-Constants.NTERM_FIX_MOD, 0, 1, PeakProperty.N_TERM_Y_ION_ONLY));

        selectedPeak = new ArrayList<Peak>(selected);       
        setScoreOfSelectedPeaks(selectedPeak, 1, Constants.massToleranceForDenovo);


        return selectedPeak.size();
    }   


//      Select Peaks only bigger than the average intensity of the given spectrum
//      input  : spectrum (already normalized)
//      output : selectedPeak size, selected peak( assign to the class value of the given spectrum )
//      2015.09.29
//      Jonghun Park

    public int selectPeakBiggerThanAvg(Spectrum spectrum) {
        // TODO Auto-generated method stub
        int i = 0;
        double sumOfIntensity = 0;
        double avgOfIntensity = 0;

        for(i = 0; i < this.size(); i++){
            sumOfIntensity += this.get(i).getNormIntensity();
        }
        avgOfIntensity = sumOfIntensity / this.size();
        //System.out.println("sum of intensity =" + sumOfIntensity);
        //System.out.println("avg of intensity =" + avgOfIntensity);

        //iterator
        //getMass()랑 그냥 Mass
        Iterator<Peak> it = this.iterator();
        Peak curPeak;
        TreeSet<Peak> selected = new TreeSet<Peak> (new MassComparator());

        while(it.hasNext()){
            curPeak = it.next();
            if(curPeak.getNormIntensity() > avgOfIntensity){
                selected.add(curPeak);
            }
        }

        selectedPeak = new ArrayList<Peak>(selected);
        return selectedPeak.size();
    }
//    /*
//     * Select Peaks only bigger than the average intensity of the peaks which have lower than 5 times of the average intensity of the given spectrum
//     * 전체 peak의 intensity 평균의 5배보다 작은 peak들의 평균보다 큰 peak들만 선택.
//     * input  : spectrum
//     * output : selectedPeak size, selected peak( assign to the class value of the given spectrum )
//     * 2015.09.29
//     * Jonghun Park
//     */
    public int selectPeakBiggerThan5Avg(Spectrum spectrum) {
        // TODO Auto-generated method stub
        int i = 0;
        double sumOfIntensity = 0;
        double avgOfIntensity = 0;

        for(i = 0; i < this.size(); i++){
            sumOfIntensity += this.get(i).getIntensity();
        }
        avgOfIntensity = sumOfIntensity / this.size();
        //System.out.println("sum of intenstiy =" + sumOfIntensity);
        //System.out.println("avg of intenstiy =" + avgOfIntensity);

        Iterator<Peak> it = this.iterator();
        Peak curPeak;
        ArrayList<Peak> firstSelected = new ArrayList<Peak>();
        TreeSet<Peak> selected = new TreeSet<Peak> (new MassComparator());

        while(it.hasNext()){
            curPeak = it.next();
            if(curPeak.getIntensity() < 5 * avgOfIntensity){
                firstSelected.add(curPeak);
            }
        }
        sumOfIntensity = 0;
        avgOfIntensity = 0;

        for(i = 0; i < firstSelected.size(); i++){
            sumOfIntensity += firstSelected.get(i).getIntensity();
        }
        avgOfIntensity = sumOfIntensity / this.size();
        //System.out.println("sum of intenstiy =" + sumOfIntensity);
        //System.out.println("avg of intenstiy =" + avgOfIntensity);

        //firstSelected mass순으로 정렬된거 보장?
        for(i = 0; i < firstSelected.size(); i++){
            if(firstSelected.get(i).getIntensity() > avgOfIntensity){
                selected.add(firstSelected.get(i));
            }
        }
        selectedPeak = new ArrayList<Peak>(selected);
        return selectedPeak.size();
    }
    public void setSelectedPeak(ArrayList<Peak> selected){
        this.selectedPeak = new ArrayList<Peak>(selected);
    }

    public ArrayList<Peak> sortbyIntensity(Spectrum spec)
    {
        TreeSet<Peak> selected = new TreeSet<Peak> (new IntensityComparator());

        for(int i = 0; i < spec.size(); i ++)
        {
            selected.add(spec.get(i));
            //System.out.println(spec.get(i).getMass());
            //System.out.println(spec.get(i).getIntensity());
        }
        ArrayList<Peak> End = new ArrayList<Peak>(selected);
        /*
        for(int i = 0; i < selected.size(); i ++)
        {
            System.out.println(End.get(i).getMass());
            System.out.println(End.get(i).getIntensity());

        }
        System.out.println(selected.size() + spec.size() + End.size());
        */
        return End;
    }
    /*
     * selected peak을 대상으로 sorting하는 method
     */
    public ArrayList<Peak> selectedSortbyIntensity(Spectrum spec)
    {
        TreeSet<Peak> selected = new TreeSet<Peak> (new IntensityComparator());
        ArrayList<Peak> preSelected = spec.getSelectedPeaks();
        for(int i = 0; i < preSelected.size() ; i ++)
        {
            selected.add(preSelected.get(i));
            //System.out.println(spec.get(i).getMass());
            //System.out.println(spec.get(i).getIntensity());
        }
        ArrayList<Peak> End = new ArrayList<Peak>(selected);
        /*
        for(int i = 0; i < selected.size(); i ++)
        {
            System.out.println(End.get(i).getMass());
            System.out.println(End.get(i).getIntensity());

        }
        System.out.println(selected.size() + spec.size() + End.size());
        */
        return End;
    }
/*    public int peakSelection (int globalSelectedSize, double selectionWindowSize, int numOfPeaksInWindow ){

        TreeSet<Peak> selected = new TreeSet<Peak> (new MassComparator());

        assert(this.checkSorted());   
        if(globalSelectedSize > this.size())   
            globalSelectedSize = this.size();
        if(this.size() == 0)   
            return -1;   

        Iterator<Peak> it = this.iterator();   
        TreeSet<Peak> globalSelected = new TreeSet<Peak>(new IntensityComparator());     

        for(int i=0; it.hasNext() && i<globalSelectedSize; i++){   
            globalSelected.add(it.next());   
        }

        TreeSet<Peak> notSelected = new TreeSet<Peak>(new MassComparator());   
        Peak curPeak;   
        while(it.hasNext())   
        {
            curPeak = it.next();   
            if(curPeak.getIntensity() > globalSelected.first().getIntensity())                                                                 
            {
                notSelected.add(globalSelected.first());   
                globalSelected.remove(globalSelected.first());   
                globalSelected.add(curPeak);

            } 
            else 
                notSelected.add(curPeak);     
        }   

        selectedPeak = new ArrayList<Peak>(globalSelected);


        return globalSelected.size();
        //return selectedPeak.size();



    }
    */
/*---------------------------------------------------------------------------------------------------------------
  ---------------------------------------------------------------------------------------------------------------
  ---------------------------------------------------------------------------------------------------------------*/
    public     void    printSelectedPeak()
    {
        System.out.println(selectedPeak.size() + " peaks are selected");
        for(Peak peak : selectedPeak)
            System.out.println(peak);
    }

    public int binarySearchforPeaks( double left ) //input : left mass
    {
        int index;   
        if( left <= this.get(0).getMass() ) //초기값 처
            index= 0;
        else if( left > this.get(this.size()-1).getMass() ) //끝처
            index= this.size()-1;
        else
        {
            int M, L= 0, R= this.size()-1;
            while( R - L > 1 )
            {
                M= ( L + R ) /2;

                if( left <= this.get(M).getMass() )
                    R= M;
                else
                    L= M;
            }
            index= R;
        }   
        return index;
    }

    public double getMatchedPeak( double mz ){
        double it=0;

        int id= binarySearchforPeaks( mz- Constants.fragmentTolerance );
        if( this.get(id).getMass() < mz - Constants.fragmentTolerance )
            return it;

        while( this.get(id).getMass() <= mz + Constants.fragmentTolerance )
        {
            if( this.get(id).getNormIntensity() > it )
                it = this.get(id).getNormIntensity();       
            id++;
            if( id == this.size() )
                break;
        }
        return it;
    }

    public double getPeakEvidence( Peak p ){
        double mz= p.getMass();   
        double ev= p.getNormIntensity(), it=0;

        if( p.property == PeakProperty.N_TERM_B_ION_ONLY || 
                p.property == PeakProperty.N_TERM_Y_ION_ONLY)
            return 1;

        if( p.property == PeakProperty.C_TERM_B_ION_ONLY || 
                p.property == PeakProperty.C_TERM_Y_ION_ONLY)
            return 1;

        if( p.getNormIntensity() < (it = getMatchedPeak(mz-1.0)) ){ return ev; }    //Isotope   

        ev += getMatchedPeak( mz+1.0 );
        if( p.getNormIntensity() >= (it=getMatchedPeak(mz-17.0)) ){ ev += it; }    //NH3
        if( p.getNormIntensity() >= (it=getMatchedPeak(mz-18.0)) ){ ev += it; }    //H2O

        ev += getMatchedPeak( p.getComplementMass(correctedMW) );         //Complementary */

        return ev;
    }

    public boolean isConfidentPeak( Peak tar, double factor ){

        int me = tar.getIndex();
        if( me < 0 ) return true;
        if( me > 0  && tar.getMass() - this.get(me-1).getMass() < 1.5 ){
            if( tar.getIntensity()/5 < this.get(me-1).getIntensity() ) return false;
        }

        if( me < this.size()-1 && this.get(me+1).getMass() - tar.getMass() < 1.5 ){
            if( tar.getIntensity() < this.get(me+1).getIntensity() ) return false;
        }           

        return true;       
    }


    public boolean checkSorted()    // identical to Tag's method
    {
        Peak tmp = null;
        for (Peak p : this)
        {
            if (tmp != null)
            {
                if (p.compareTo(tmp)<0) return false;
            }
            tmp = p;
        }
        return true;
    }

    public int compareTo(Spectrum s) 
    {
        if( observedMW > s.observedMW ) return 1;
        else if( observedMW == s.observedMW ) return 0;
        else return -1;
    }


    private void setScoreOfSelectedPeaks( ArrayList<Peak> selected, int assumedCS, double tolerance )
    {
        double isotopeDelta = 1.0034, NH3Delta = Constants.NH3, H2ODelta = Constants.H2O;
        for( int node=0; node<selected.size(); node++){

            if( selected.get(node).property == PeakProperty.N_TERM_B_ION_ONLY || 
                    selected.get(node).property == PeakProperty.N_TERM_Y_ION_ONLY ){
                selected.get(node).setProbability(1.);
                continue;
            }

            if( selected.get(node).property == PeakProperty.C_TERM_B_ION_ONLY || 
                    selected.get(node).property == PeakProperty.C_TERM_Y_ION_ONLY ){
                selected.get(node).setProbability(1.);
                continue;
            }

            double targetmz=0;       
            double score = selected.get(node).getNormIntensity(); 
            double ionMZ= selected.get(node).getMass(); 
            double ionIT= selected.get(node).getIntensity();

            // check this is isotope???
            int OK = -1;       
            double okmax = 0;
            targetmz= ionMZ-isotopeDelta;
            for(int i=node-1; i>-1; i--){
                if( selected.get(i).getMass() < targetmz-tolerance ) break;           
                else if( Math.abs(selected.get(i).getMass()-targetmz) < tolerance ){
                    if( selected.get(i).getIntensity() > okmax ) {
                        okmax = selected.get(i).getIntensity();
                        OK = i;
                    }
                }
            }
            if( OK != -1 ) {
                if( okmax > ionIT ) selected.get(node).setProbability(0.);
                else selected.get(node).setProbability(score);
                continue;
            }//*/


            int ISO = node + 1; // plus isotope peak   
            double prevISO = ionIT;
            boolean isoDecent = false;
            targetmz= ionMZ + isotopeDelta;
            for( int nth_iso=1 ;  ; nth_iso++  ){
                int H = -1;
                double imax = 0;               
                for(int i=ISO; i<selected.size(); i++){
                    if( selected.get(i).getMass() > targetmz+Constants.massToleranceForDenovo ) break;           
                    else if( Math.abs(selected.get(i).getMass()-targetmz) < Constants.massToleranceForDenovo ){
                        if( selected.get(i).getIntensity() > imax ) {
                            imax = selected.get(i).getIntensity();
                            H = i;
                        }
                    }
                }

            //    if( H == -1 ) break;
                if( H == -1 || prevISO < imax/5 ) break;
                if( isoDecent && prevISO < imax ) break;
                if( prevISO > imax ) isoDecent=true;

                score += selected.get(H).getNormIntensity()/nth_iso;   
                ISO = H+1;
                targetmz= selected.get(H).getMass() + isotopeDelta;       
                prevISO = imax;
            }//*/

            int NLOSS = -1;       
            double lossmax=0, lossmz=0;
            targetmz= ionMZ-H2ODelta;
            for(int i=node-1; i>-1; i--){
                if( selected.get(i).getMass() < targetmz-tolerance ) break;           
                else if( Math.abs(selected.get(i).getMass()-(ionMZ-NH3Delta)) < tolerance ){
                    if( selected.get(i).getIntensity() > lossmax ) {
                        lossmax = selected.get(i).getIntensity();
                        lossmz = selected.get(i).getMass();
                        NLOSS = i;
                    }
                }
                else if( Math.abs(selected.get(i).getMass()-targetmz) < tolerance ){
                    if( selected.get(i).getIntensity() > lossmax ) {
                        lossmax = selected.get(i).getIntensity();
                        lossmz = selected.get(i).getMass();
                        NLOSS = i;
                    }
                }
            }
            if( NLOSS != -1 ){
                double lossScore = selected.get(NLOSS).getNormIntensity();               
                int NL_H = -1; // plus isotope peak       
                double nsmax = 0;
                targetmz= lossmz+isotopeDelta;
                for(int i=NLOSS+1; i<selected.size(); i++){
                    if( selected.get(i).getMass() > targetmz+tolerance ) break;           
                    else if( Math.abs(selected.get(i).getMass()-targetmz) < tolerance ){
                        if( selected.get(i).getIntensity() > nsmax ) {
                            nsmax = selected.get(i).getIntensity();
                            NL_H = i;
                        }
                    }
                }
                if( NL_H != -1 ) lossScore += selected.get(NL_H).getNormIntensity();               
                score += lossScore*0.5;

            }//*/

            double complement = correctedMW - ionMZ + Constants.Proton*2;   
            score += getMatchedPeak( complement ); 

            selected.get(node).setProbability(score);
        }
    }

}

class PrecursorComparator implements Comparator<Spectrum>
{
    public int compare(Spectrum spec1, Spectrum spec2)
    {
        double preMass1 = spec1.getPrecursor();
        double preMass2 = spec2.getPrecursor();
        if(preMass1 > preMass2)
            return 1;
        else if(preMass1 == preMass2)
            return 0;
        else 
            return -1;
    }

    public boolean equals(Spectrum spec1, Spectrum spec2)
    {
        return spec1.getPrecursor() == spec2.getPrecursor();
    }
}

class SpecNameComparator implements Comparator<Spectrum>
{
    public int compare(Spectrum spec1, Spectrum spec2)
    {   
        if( spec1.getName().compareTo(spec2.getName()) > 0 )
            return 1;
        else if( spec1.getName().compareTo(spec2.getName()) == 0 )
            return 0;
        else 
            return -1;
    }

    public boolean equals(Spectrum spec1, Spectrum spec2)
    {
        return spec1.getName().compareTo(spec2.getName()) == 0;
    }
}
