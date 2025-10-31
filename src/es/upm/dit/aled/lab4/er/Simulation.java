package es.upm.dit.aled.lab4.er;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Scanner;

import es.upm.dit.aled.lab4.er.gui.EmergencyRoomGUI;

/**
 * Starts a simulation of a hospital ER modeled by the EmergencyRoom class. Uses
 * EmergencyRoomGUI as a graphical user interface.
 * 
 * @author rgarciacarmona
 */
public class Simulation {

	private EmergencyRoom er;

	public static void main(String[] args) {
		try {
			Simulation sim = new Simulation(args[0]);
			sim.start();
		} catch (FileNotFoundException e) {
			System.out.println("File not found.");
		}
		// El primer argumento debe ser el nombre del archivo de datos del 
		// servicio de urgencias, es decir: er_data.txt.
	}

	/**
	 * Simulates an emergency room from a provided file
	 */
	public Simulation(String file) throws FileNotFoundException {
		// Initializes EmergencyRoom
		this.er = new EmergencyRoom();
		// Loads data from file
		File inputFile = new File(file);
		er.readFile(new Scanner(inputFile));
		// Initializes the GUI
		EmergencyRoomGUI.initialize(er);
	}

	/**
	 * Starts all Patients and waits until all of them have been discharged.
	 */
	public void start() {
		System.out.println("Starting simulation... Admitting " + er.getPatients().size() + " patients.");
		// Makes a copy of the patients and shuffle them to impose
		// a random order of arrival
		List<Patient> patientsList = new ArrayList<Patient>(er.getPatients());
		Collections.shuffle(patientsList);
		// Admit all patients
		for (Patient p : patientsList) {
			er.admit(p);
		}
		// Wait until all patient threads have finished
		for (Patient p : patientsList) {
			er.waitForDischarge(p);
		}
		System.out.println("All patients processed. Simulation finished.");
	}
}