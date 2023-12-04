import curses
from Bio import PDB
import random
import math

# Global variables for translation and mouse selection
translate_x = 0
translate_y = 0
selected_atom_info = None  # Store information about the selected atom

def load_covid_structure():
    pdb_file_path = "7p35-pdb-bundle1.pdb"
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('COVID-19', pdb_file_path)
    return structure

def generate_atom_colors(structure):
    atom_color_map = {}
    color_count = 1
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_name = atom.get_name()
                    if atom_name not in atom_color_map:
                        atom_color_map[atom_name] = color_count
                        color_count += 1
    return atom_color_map

def get_atom_description(atom_name):
    atom_descriptions = {
        'C': 'Carbon - Essential building block for organic molecules. It forms the backbone of many biological molecules.',
        'N': 'Nitrogen - Commonly found in amino groups of amino acids. Essential for the structure of proteins and nucleic acids.',
        'O': 'Oxygen - Important for cellular respiration and energy production. Found in water and many organic molecules.',
        'H': 'Hydrogen - The simplest and most abundant element in the universe. Commonly found in water and organic compounds.',
        'S': 'Sulfur - Essential for the structure of some amino acids and vitamins. Involved in disulfide bonds in proteins.',
        'P': 'Phosphorus - Key component of nucleic acids (DNA and RNA) and ATP, the primary energy currency of cells.',
        'NA': 'Sodium - Important for maintaining osmotic balance and electrical excitability in cells.',
        'CL': 'Chlorine - Essential for maintaining osmotic balance and is a component of the digestive enzyme pepsin.',
        'K': 'Potassium - Critical for nerve and muscle function. Helps maintain the body\'s fluid balance.',
        'MG': 'Magnesium - Essential for the function of many enzymes and the stability of nucleic acids and proteins.',
        'CA': 'Calcium - Important for the structure of bones and teeth. Essential for blood clotting and muscle function.',
        'FE': 'Iron - Involved in oxygen transport in hemoglobin and electron transport in cellular respiration.',
        'ZN': 'Zinc - Important for the function of many enzymes and essential for the immune system.',
        'CU': 'Copper - Necessary for the function of certain enzymes, including those involved in iron metabolism.',
        'MN': 'Manganese - Required for the activity of some enzymes involved in metabolism and bone formation.',
        'I': 'Iodine - Essential for the synthesis of thyroid hormones, which regulate metabolism.',
        'F': 'Fluorine - Strengthens tooth enamel and is used in some dental products.',
        'CB': 'Beta Carbon - Found in amino acids. Involved in the side chain of amino acids.',
        'CG': 'Gamma Carbon - Found in amino acids. Involved in the side chain of amino acids.',
        'CD': 'Delta Carbon - Found in amino acids. Involved in the side chain of amino acids.',
        'CE': 'Epsilon Carbon - Found in amino acids. Involved in the side chain of amino acids.',
        'NZ': 'Epsilon Nitrogen - Found in amino acids. Part of the side chain of lysine.',
        'A': 'Adenine - One of the four nucleobases in DNA and RNA. Forms base pairs with thymine (in DNA) or uracil (in RNA).',
        'C5M': '5-Methylcytidine - Modified nucleoside found in RNA.',
        'G': 'Guanine - One of the four nucleobases in DNA and RNA. Forms base pairs with cytosine.',
        'U': 'Uracil - One of the four nucleobases in RNA. Replaces thymine in the RNA structure.',
        'RNA': 'Ribonucleic Acid - Genetic material of many viruses, including SARS-CoV-2.',
        'CO3': 'Carbonate - Component of some biological molecules and buffers.',
        'PO4': 'Phosphate - Essential for the structure of nucleic acids (DNA and RNA) and ATP.',
        'RU5P': 'Ribulose-5-Phosphate - An intermediate in the pentose phosphate pathway.',
        'PC': 'Phosphatidylcholine - A major component of cell membranes.',
        # Add more atom descriptions as needed
    }
    return atom_descriptions.get(atom_name, f"No description available for {atom_name}.")


def draw_line(win, y1, x1, y2, x2, char, color_pair):
    # Draw a line between two points on the screen
    dx = x2 - x1
    dy = y2 - y1
    steps = max(abs(dx), abs(dy)) or 1  # Ensure steps is at least 1

    for i in range(steps + 1):
        x = int(x1 + i * dx / steps)
        y = int(y1 + i * dy / steps)
        try:
            win.addch(y, x, char, color_pair)
        except curses.error:
            pass  # Ignore errors when trying to draw outside the window boundaries

def draw_bordered_text(win, y, x, text, attributes, color_pair):
    # Draw a bordered text on the screen
    border_char = '*'
    win.attron(color_pair)

    # Draw the top border
    win.addstr(y, x, border_char * (len(text) + 2), attributes)

    # Draw the text with side borders
    win.addstr(y + 1, x, f"{border_char}{text}{border_char}", attributes)

    # Draw the bottom border
    win.addstr(y + 2, x, border_char * (len(text) + 2), attributes)

    win.attroff(color_pair)

def center_on_atom(atom):
    global translate_x, translate_y

    # Set translation values to center the screen on the selected atom
    translate_x = -atom.get_coord()[0]
    translate_y = -atom.get_coord()[1]

def get_adjacent_atom(current_atom, direction, structure):
    # Get the adjacent atom based on the current focus and the direction
    for model in structure:
        for chain in model:
            for residue in chain:
                atoms = list(residue)
                try:
                    index = atoms.index(current_atom)
                    if direction == "left" and index > 0:
                        return atoms[index - 1]
                    elif direction == "right" and index < len(atoms) - 1:
                        return atoms[index + 1]
                except ValueError:
                    pass
    return current_atom

def random_atom_info(structure):
    # Select a random atom and return its information
    random_model = random.choice(structure)
    random_chain = random.choice(list(random_model))
    random_residue = random.choice(list(random_chain))
    random_atom = random.choice(list(random_residue))
    return random_atom

def draw_covid_19(structure, win, atom_color_map):
    center_x = curses.COLS // 2
    center_y = curses.LINES // 2

    for model in structure:
        for chain in model:
            for residue in chain:
                atoms = list(residue)
                for i in range(len(atoms) - 1):
                    x1, y1, _ = atoms[i].get_coord()
                    x2, y2, _ = atoms[i + 1].get_coord()

                    # Adjust coordinates to center the screen
                    x1 += center_x + translate_x
                    y1 += center_y + translate_y
                    x2 += center_x + translate_x
                    y2 += center_y + translate_y

                    # Draw lines between atoms
                    try:
                        draw_line(win, int(y1), int(x1), int(y2), int(x2), ' ', curses.color_pair(atom_color_map[atoms[i].get_name()]))
                    except curses.error:
                        pass  # Ignore errors when trying to draw outside the window boundaries

                for atom in atoms:
                    x, y, _ = atom.get_coord()

                    # Adjust coordinates to center the screen
                    x += center_x + translate_x
                    y += center_y + translate_y

                    # Draw only the first letter of the atom name
                    try:
                        win.addch(int(y), int(x), atom.get_name()[0], curses.color_pair(atom_color_map[atom.get_name()]))
                    except curses.error:
                        pass  # Ignore errors when trying to draw outside the window boundaries

    # Draw box in the lower-left corner with the full atom description
    box_height = 5
    box_width = curses.COLS // 3
    box_y = curses.LINES - box_height
    box_x = 0
    win.box()
    
    # Get the full description for the selected atom
    selected_atom_description = get_atom_description(selected_atom_info.get_name())
    
    # Draw the full description inside the box
    for i, line in enumerate(selected_atom_description.split('\n')):
        win.addstr(box_y + i + 1, box_x + 1, line, curses.color_pair(atom_color_map[selected_atom_info.get_name()]))

def main(stdscr):
    global translate_x, translate_y, selected_atom_info

    curses.curs_set(0)  # Hide the cursor
    stdscr.clear()
    stdscr.refresh()

    stdscr.addstr(0, 0, "Press 'q' to quit | Use arrow keys to navigate", curses.A_BOLD)

    covid_structure = load_covid_structure()
    atom_color_map = generate_atom_colors(covid_structure)

    curses.start_color()
    for i in range(1, len(atom_color_map) + 1):
        # Use color pairs with better visibility
        curses.init_pair(i, i % 7 + 1, curses.COLOR_BLACK)

    selected_atom_info = random_atom_info(covid_structure)
    center_on_atom(selected_atom_info)

    while True:
        key = stdscr.getch()

        if key == ord('q'):
            break
        elif key == curses.KEY_LEFT:
            selected_atom_info = get_adjacent_atom(selected_atom_info, "left", covid_structure)
            center_on_atom(selected_atom_info)
        elif key == curses.KEY_RIGHT:
            selected_atom_info = get_adjacent_atom(selected_atom_info, "right", covid_structure)
            center_on_atom(selected_atom_info)
        elif key == curses.KEY_UP:
            selected_atom_info = get_adjacent_atom(selected_atom_info, "up", covid_structure)
            center_on_atom(selected_atom_info)
        elif key == curses.KEY_DOWN:
            selected_atom_info = get_adjacent_atom(selected_atom_info, "down", covid_structure)
            center_on_atom(selected_atom_info)

        stdscr.clear()

        # Draw COVID-19 structure
        draw_covid_19(covid_structure, stdscr, atom_color_map)

        # Draw selected atom info with a rectangle, larger text, and a border
        if selected_atom_info:
            atom_info = f"Atom: {selected_atom_info.get_name()} in {selected_atom_info.parent.resname} (Chain {selected_atom_info.parent.parent.id})"
            draw_bordered_text(stdscr, curses.LINES - 3, 0, atom_info, curses.A_BOLD, curses.color_pair(atom_color_map[selected_atom_info.get_name()]))

        stdscr.refresh()

if __name__ == "__main__":
    curses.wrapper(main)
