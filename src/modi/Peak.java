package modi;

import java.util.Comparator;

public class Peak implements Comparable<Peak> {
  double mass;
  double intensity;
  double normIntensity = 0;
  int charge = 1;
  int index = -1;
  final PeakProperty property;

  public Peak(int index, double mass, double intensity, int charge, PeakProperty property) {
    this.index = index;
    this.mass = mass;
    this.intensity = intensity;
    this.charge = charge;
    this.property = property;
  }

  public Peak(int index, double mass, double intensity, int charge) {
    this(index, mass, intensity, charge, PeakProperty.NORMAL);
  }

  public Peak(int index, double mass, double intensity) {
    this(index, mass, intensity, 1, PeakProperty.NORMAL);
  }

  public int getIndex() {
    return index;
  }

  public double getMass() {
    return mass;
  }

  public double getIntensity() {
    return intensity;
  }

  public double getNormIntensity() {
    return normIntensity;
  }

  public int getCharge() {
    return charge;
  }

  public PeakProperty getPeakProperty() {
    return property;
  }

  public double getComplementMass(double motherMass) {
    return motherMass - mass + 2 * Constants.Proton;
  }

  public void set(double m, double i) {
    mass = m;
    intensity = i;
  }

  public void shiftMass(double m) {
    mass += m;
  }

  // normalized intensity, and charge may be decided after Peak is made
  public void setNormIntensity(double normIntensity) {
    this.normIntensity = normIntensity;
  }

  public void setCharge(int charge) {
    this.charge = charge;
  }

  public Peak getShiftedPeak(double shiftMass) {
    Peak p = this.clone();
    p.mass += shiftMass;

    return p;
  }

  public int compareTo(Peak p) {   // default comparator : mass
    if (mass > p.mass)
      return 1;
    else if (mass < p.mass)
      return -1;
    else
      return 0;
  }

  public static double getMassDifference(Peak p1, Peak p2) {
    return Math.abs(p1.mass - p2.mass);
  }

  public String toString() {
    return new String(Constants.getString(mass) + " " + this.property);
  }

  public String toString2() {
    return Constants.getString(mass) + " "
         + Constants.getString(normIntensity) + '\n';
  }

  public Peak clone() {
    Peak p = new Peak(index, mass, intensity, charge, property);
    p.setNormIntensity(normIntensity);

    return p;
  }

  private double probability = 0;

  public double getProbability() {
    return probability;
  }

  public void setProbability(double pa) {
    probability = pa;
  }
}

class MassComparator implements Comparator<Peak> {
  public int compare(Peak p1, Peak p2) {
    if (p1.mass > p2.mass)
      return 1;
    else if (p1.mass == p2.mass) {
      if (p1.intensity > p2.intensity)
        return 1;
      else if (p1.intensity == p2.intensity)
        return 0;
      else
        return -1;
    }
    else
      return -1;
  }

  public boolean equals(Peak p1, Peak p2) {
    return p1.mass == p2.mass && p1.intensity == p2.intensity;
  }
}

class IntensityComparator implements Comparator<Peak> {
  public int compare(Peak p1, Peak p2) {
    if (p1.intensity > p2.intensity)
      return 1;
    else if (p1.intensity == p2.intensity) {
      if (p1.mass > p2.mass)
        return 1;
      else if (p1.mass == p2.mass)
        return 0;
      else
        return -1;
    }
    else
      return -1;  //*/
  }

  public boolean equals(Peak p1, Peak p2) {
    return p1.mass == p2.mass && p1.intensity == p2.intensity;
  }
}