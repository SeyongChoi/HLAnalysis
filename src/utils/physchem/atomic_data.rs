// atomic number(int) -> input: symbol(String)
// atomic mass(f64) -> input: symobl(String) or atomic number(int)
use pyo3::prelude::*;
use std::collections::HashMap;

// Atomic data
lazy_static::lazy_static! {
    static ref ELEMENTS: HashMap<String, (i32, f64)> = {
        let mut m = HashMap::new();
        m.insert("H".to_string(), (1, 1.008));    // Hydrogen
        m.insert("He".to_string(), (2, 4.0026));  // Helium
        m.insert("Li".to_string(), (3, 6.94));    // Lithium
        m.insert("Be".to_string(), (4, 9.0122));  // Beryllium
        m.insert("B".to_string(), (5, 10.81));    // Boron
        m.insert("C".to_string(), (6, 12.011));   // Carbon
        m.insert("N".to_string(), (7, 14.007));   // Nitrogen
        m.insert("O".to_string(), (8, 15.999));   // Oxygen
        m.insert("F".to_string(), (9, 18.998));   // Fluorine
        m.insert("Ne".to_string(), (10, 20.180)); // Neon
        m.insert("Na".to_string(), (11, 22.990)); // Sodium
        m.insert("Mg".to_string(), (12, 24.305)); // Magnesium
        m.insert("Al".to_string(), (13, 26.982)); // Aluminium
        m.insert("Si".to_string(), (14, 28.085)); // Silicon
        m.insert("P".to_string(), (15, 30.974));  // Phosphorus
        m.insert("S".to_string(), (16, 32.06));   // Sulfur
        m.insert("Cl".to_string(), (17, 35.45));  // Chlorine
        m.insert("Ar".to_string(), (18, 39.948)); // Argon
        m.insert("K".to_string(), (19, 39.098));  // Potassium
        m.insert("Ca".to_string(), (20, 40.078)); // Calcium
        m.insert("Sc".to_string(), (21, 44.956)); // Scandium
        m.insert("Ti".to_string(), (22, 47.867)); // Titanium
        m.insert("V".to_string(), (23, 50.942));  // Vanadium
        m.insert("Cr".to_string(), (24, 52.0));   // Chromium
        m.insert("Mn".to_string(), (25, 54.938)); // Manganese
        m.insert("Fe".to_string(), (26, 55.845)); // Iron
        m.insert("Ni".to_string(), (27, 58.693)); // Nickel
        m.insert("Co".to_string(), (28, 58.933)); // Cobalt
        m.insert("Cu".to_string(), (29, 63.546)); // Copper
        m.insert("Zn".to_string(), (30, 65.38));  // Zinc
        m.insert("Ga".to_string(), (31, 69.723)); // Gallium
        m.insert("Ge".to_string(), (32, 72.63));  // Germanium
        m.insert("As".to_string(), (33, 74.922)); // Arsenic
        m.insert("Se".to_string(), (34, 78.971)); // Selenium
        m.insert("Br".to_string(), (35, 79.904)); // Bromine
        m.insert("Kr".to_string(), (36, 83.798)); // Krypton
        m.insert("Rb".to_string(), (37, 85.468)); // Rubidium
        m.insert("Sr".to_string(), (38, 87.62));  // Strontium
        m.insert("Y".to_string(), (39, 88.906));  // Yttrium
        m.insert("Zr".to_string(), (40, 91.224)); // Zirconium
        m.insert("Nb".to_string(), (41, 92.906)); // Niobium
        m.insert("Mo".to_string(), (42, 95.95));  // Molybdenum
        m.insert("Tc".to_string(), (43, 98.0));   // Technetium
        m.insert("Ru".to_string(), (44, 101.07)); // Ruthenium
        m.insert("Rh".to_string(), (45, 102.91)); // Rhodium
        m.insert("Pd".to_string(), (46, 106.42)); // Palladium
        m.insert("Ag".to_string(), (47, 107.87)); // Silver
        m.insert("Cd".to_string(), (48, 112.41)); // Cadmium
        m.insert("In".to_string(), (49, 114.82)); // Indium
        m.insert("Sn".to_string(), (50, 118.71)); // Tin
        m.insert("Sb".to_string(), (51, 121.76)); // Antimony
        m.insert("I".to_string(), (53, 126.90));  // Iodine
        m.insert("Te".to_string(), (52, 127.60)); // Tellurium
        m.insert("Xe".to_string(), (54, 131.29)); // Xenon
        m.insert("Cs".to_string(), (55, 132.91)); // Cesium
        m.insert("Ba".to_string(), (56, 137.33)); // Barium
        m.insert("La".to_string(), (57, 138.91)); // Lanthanum
        m.insert("Ce".to_string(), (58, 140.12)); // Cerium
        m.insert("Pr".to_string(), (59, 140.91)); // Praseodymium
        m.insert("Nd".to_string(), (60, 144.24)); // Neodymium
        m.insert("Pm".to_string(), (61, 145.0));  // Promethium
        m.insert("Sm".to_string(), (62, 150.36)); // Samarium
        m.insert("Eu".to_string(), (63, 152.00)); // Europium
        m.insert("Gd".to_string(), (64, 157.25)); // Gadolinium
        m.insert("Tb".to_string(), (65, 158.93)); // Terbium
        m.insert("Dy".to_string(), (66, 162.50)); // Dysprosium
        m.insert("Ho".to_string(), (67, 164.93)); // Holmium
        m.insert("Er".to_string(), (68, 167.26)); // Erbium
        m.insert("Tm".to_string(), (69, 168.93)); // Thulium
        m.insert("Yb".to_string(), (70, 173.04)); // Ytterbium
        m.insert("Lu".to_string(), (71, 175.00)); // Lutetium
        m.insert("Hf".to_string(), (72, 178.49)); // Hafnium
        m.insert("Ta".to_string(), (73, 180.95)); // Tantalum
        m.insert("W".to_string(), (74, 183.84));  // Tungsten
        m.insert("Re".to_string(), (75, 186.21)); // Rhenium
        m.insert("Os".to_string(), (76, 190.23)); // Osmium
        m.insert("Ir".to_string(), (77, 192.22)); // Iridium
        m.insert("Pt".to_string(), (78, 195.08)); // Platinum
        m.insert("Au".to_string(), (79, 196.97)); // Gold
        m.insert("Hg".to_string(), (80, 200.59)); // Mercury
        m.insert("Tl".to_string(), (81, 204.38)); // Thallium
        m.insert("Pb".to_string(), (82, 207.2));  // Lead
        m.insert("Bi".to_string(), (83, 208.98)); // Bismuth
        m.insert("Po".to_string(), (84, 209.0));  // Polonium
        m.insert("At".to_string(), (85, 210.0));  // Astatine
        m.insert("Rn".to_string(), (86, 222.0));  // Radon
        m.insert("Fr".to_string(), (87, 223.0));  // Francium
        m.insert("Ra".to_string(), (88, 226.0));  // Radium
        m.insert("Ac".to_string(), (89, 227.0));  // Actinium
        m.insert("Th".to_string(), (90, 232.04)); // Thorium
        m.insert("Pa".to_string(), (91, 231.04)); // Protactinium
        m.insert("U".to_string(), (92, 238.03));  // Uranium
        m.insert("Np".to_string(), (93, 237.0));  // Neptunium
        m.insert("Pu".to_string(), (94, 244.0));  // Plutonium
        m.insert("Am".to_string(), (95, 243.0));  // Americium
        m.insert("Cm".to_string(), (96, 247.0));  // Curium
        m.insert("Bk".to_string(), (97, 247.0));  // Berkelium
        m.insert("Cf".to_string(), (98, 251.0));  // Californium
        m.insert("Es".to_string(), (99, 252.0));  // Einsteinium
        m.insert("Fm".to_string(), (100, 257.0)); // Fermium
        m.insert("Md".to_string(), (101, 258.0)); // Mendelevium
        m.insert("No".to_string(), (102, 259.0)); // Nobelium
        m.insert("Lr".to_string(), (103, 262.0)); // Lawrencium
        m.insert("Rf".to_string(), (104, 267.0)); // Rutherfordium
        m.insert("Db".to_string(), (105, 270.0)); // Dubnium
        m.insert("Sg".to_string(), (106, 271.0)); // Seaborgium
        m.insert("Bh".to_string(), (107, 270.0)); // Bohrium
        m.insert("Hs".to_string(), (108, 277.0)); // Hassium
        m.insert("Mt".to_string(), (109, 276.0)); // Meitnerium
        m.insert("Ds".to_string(), (110, 281.0)); // Darmstadtium
        m.insert("Rg".to_string(), (111, 280.0)); // Roentgenium
        m.insert("Cn".to_string(), (112, 285.0)); // Copernicium
        m.insert("Nh".to_string(), (113, 284.0)); // Nihonium
        m.insert("Fl".to_string(), (114, 289.0)); // Flerovium
        m.insert("Mc".to_string(), (115, 288.0)); // Moscovium
        m.insert("Lv".to_string(), (116, 293.0)); // Livermorium
        m.insert("Ts".to_string(), (117, 294.0)); // Tennessine
        m.insert("Og".to_string(), (118, 294.0)); // Oganesson
        m
    };
}

// Function to get atomic number from symbol
#[pyfunction]
pub fn atomic_numbers(symbol: String) -> PyResult<i32> {
    match ELEMENTS.get(&symbol) {
        Some((number, _)) => Ok(*number),
        None => Err(pyo3::exceptions::PyKeyError::new_err("Element not found")),
    }
}
// Function to get atomic mass from symbol or atomic number
#[pyfunction]
pub fn atomic_masses(py: Python, input: PyObject) -> PyResult<f64> {
    // Try to extract the input as a String
    if let Ok(symbol) = input.extract::<String>(py) {
        // Match the symbol in the ELEMENTS map
        match ELEMENTS.get(&symbol) {
            Some((_, mass)) => Ok(*mass),
            None => Err(pyo3::exceptions::PyKeyError::new_err("Element not found")),
        }
    } else if let Ok(number) = input.extract::<i32>(py) {
        // If input is atomic number (i32)
        for (_, (num, mass)) in ELEMENTS.iter() {
            if *num == number {
                return Ok(*mass);
            }
        }
        Err(pyo3::exceptions::PyKeyError::new_err("Element not found"))
    } else {
        // If the input is neither a symbol (String) nor an atomic number (i32)
        Err(pyo3::exceptions::PyTypeError::new_err("Invalid input type"))
    }
}
