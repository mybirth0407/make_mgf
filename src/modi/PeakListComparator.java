package modi;

import java.util.ArrayList;
import java.util.Iterator;

// ArrayList<Peak> must be sorted by peak mass
public class PeakListComparator {
  private ArrayList<PeakPair> pairedPeaks = new ArrayList<PeakPair>();
  private ArrayList<Peak> firstOnly = new ArrayList<Peak>();
  private ArrayList<Peak> secondOnly = new ArrayList<Peak>();
  private ArrayList<PeakPair> shared = new ArrayList<PeakPair>();

  public PeakListComparator(ArrayList<Peak> firstPeakList, ArrayList<Peak> secondPeakList) {
    Iterator<Peak> it1 = firstPeakList.iterator();
    Iterator<Peak> it2 = secondPeakList.iterator();
    Peak p1 = null, p2 = null;
    boolean p1Inserted = true, p2Inserted = true;
    while (it1.hasNext() || it2.hasNext()) {
      if (p1Inserted) {
        if (!it1.hasNext())
          break;
        p1 = it1.next();
        p1Inserted = false;
      }
      if (p2Inserted) {
        if (!it2.hasNext())
          break;
        p2 = it2.next();
        p2Inserted = false;
      }

      if (Constants.fEqual(p1.getMass(), p2.getMass())) {
        pairedPeaks.add(new PeakPair(p1, p2));
        shared.add(new PeakPair(p1, p2));
        p1Inserted = p2Inserted = true;
      }
      else if (p1.getMass() < p2.getMass() - Constants.fragmentTolerance) {
        pairedPeaks.add(new PeakPair(p1, null));
        firstOnly.add(p1);
        p1Inserted = true;
      }
      else {
        pairedPeaks.add(new PeakPair(null, p2));
        secondOnly.add(p2);
        p2Inserted = true;
      }
    }
    if (p1Inserted == false) {
      firstOnly.add(p1);
      pairedPeaks.add(new PeakPair(p1, null));
    }
    if (p2Inserted == false) {
      secondOnly.add(p2);
      pairedPeaks.add(new PeakPair(null, p2));
    }

    while (it1.hasNext()) {
      p1 = it1.next();
      firstOnly.add(p1);
      pairedPeaks.add(new PeakPair(p1, null));
    }
    while (it2.hasNext()) {
      p2 = it2.next();
      secondOnly.add(p2);
      pairedPeaks.add(new PeakPair(null, p2));
    }
  }

  public ArrayList<Peak> getFirstOnly() {
    return firstOnly;
  }

  public ArrayList<Peak> getSecondOnly() {
    return secondOnly;
  }

  public ArrayList<PeakPair> getPairedPeaks() {
    return pairedPeaks;
  }

  public ArrayList<PeakPair> getSharedPeaks() {
    return shared;
  }

  public int getSharedPeakCount() {
    return shared.size();
  }

  public double getSharedPeakIntensity(int mode) {
    // mode 0 : intensity sum
    // mode 1 : more intense peak sum
    // mode 2 : less intense peak sum
    double intensitySum = 0.;
    for (PeakPair pair : shared) {
      double first = pair.getFirst().getNormIntensity();
      double second = pair.getSecond().getNormIntensity();
      if (mode == 0)
        intensitySum += first + second;
      else if (mode == 1) {
        intensitySum += first > second ? first : second;
      }
      else if (mode == 2) {
        intensitySum += first > second ? second : first;
      }
    }
    return intensitySum;
  }

  public ArrayList<Peak> getMergedPeaks()    // consist of high intensity peaks
  {
    ArrayList<Peak> mergedPeaks = new ArrayList<Peak>();
    ArrayList<PeakPair> pairedList = getPairedPeaks();
    Iterator<PeakPair> it = pairedList.iterator();
    while (it.hasNext()) {
      PeakPair peakPair = it.next();
      if (peakPair.first == null)
        mergedPeaks.add(peakPair.second);
      else if (peakPair.second == null)
        mergedPeaks.add(peakPair.first);
      else {
        if (peakPair.first.getIntensity() > peakPair.second.getIntensity())
          mergedPeaks.add(peakPair.first);
        else
          mergedPeaks.add(peakPair.second);
      }
    }

    return mergedPeaks;
  }

  public boolean firstContainsSecond() {
    return secondOnly.size() == 0;
  }

  public boolean secondContainsFirst() {
    return firstOnly.size() == 0;
  }

  // merge second peak list to first peak list
  public static void mergePeakList(ArrayList<Peak> firstPeakList, ArrayList<Peak> secondPeakList) {
    ArrayList<Peak>
        mergedPeakList =
        new PeakListComparator(firstPeakList, secondPeakList).getMergedPeaks();
    firstPeakList.clear();
    firstPeakList.addAll(mergedPeakList);
  }
}