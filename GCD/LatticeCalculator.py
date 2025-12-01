import numpy as np


def generateLatticeSpots(nstrings=10):
    offset = np.pi * (1.0 / 6)
    anglediff = np.pi * (1.0 / 3)
    neighbourangles = [
        offset,
        anglediff + offset,
        2.0 * anglediff + offset,
        3.0 * anglediff + offset,
        4.0 * anglediff + offset,
        5.0 * anglediff + offset,
    ]

    stringposx = [0.0]
    stringposy = [0.0]
    theta = [0.0]
    spacing = 50.0

    while len(stringposx) < nstrings:
        minradius = 1000000.0
        minradstring = 0
        minradstringneighbours = 10000
        for i in range(len(stringposx)):
            nneighbours = 0
            rad = np.sqrt((stringposx[i]) ** 2.0 + (stringposy[i]) ** 2.0)
            for j in range(len(stringposx)):
                if i == j:
                    continue
                dist = np.sqrt(
                    (stringposx[j] - stringposx[i]) ** 2.0
                    + (stringposy[j] - stringposy[i]) ** 2.0
                )
                if dist < spacing * 1.2:
                    nneighbours += 1
            if nneighbours < len(neighbourangles) and rad <= minradius:
                if rad < minradius:
                    minradius = rad
                    minradstring = i
                    minradstringneighbours = nneighbours
                elif nneighbours < minradstringneighbours:
                    minradius = rad
                    minradstring = i
                    minradstringneighbours = nneighbours

        maxneighours = 0
        maxneighbourstring = 0
        for j in range(len(neighbourangles)):
            newposx = stringposx[minradstring] + spacing * np.sin(neighbourangles[j])
            newposy = stringposy[minradstring] + spacing * np.cos(neighbourangles[j])

            nneighbours = 0
            overlap = False
            for k in range(len(stringposx)):
                dist = np.sqrt(
                    (newposx - stringposx[k]) ** 2.0 + (newposy - stringposy[k]) ** 2.0
                )
                if dist < spacing * 0.8:
                    nneighbours = 0
                    overlap = True
                if dist < spacing * 1.2:
                    nneighbours += 1
            if nneighbours > maxneighours and not overlap:
                maxneighours = nneighbours
                maxneighbourstring = j
        stringposx.append(
            stringposx[minradstring]
            + spacing * np.sin(neighbourangles[maxneighbourstring])
        )
        stringposy.append(
            stringposy[minradstring]
            + spacing * np.cos(neighbourangles[maxneighbourstring])
        )
        theta.append(neighbourangles[maxneighbourstring])

    mean_x = sum(stringposx) / len(stringposx)
    mean_y = sum(stringposy) / len(stringposy)

    for i in range(len(stringposx)):
        stringposx[i] = (stringposx[i] - mean_x) / spacing
        stringposy[i] = (stringposy[i] - mean_y) / spacing

    return stringposx, stringposy, theta

#Generated main function to show usage of generateLatticeSpots
if __name__ == "__main__":
    """Main function to demonstrate usage of generateLatticeSpots."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate lattice structure for string positions")
    parser.add_argument("-n", "--nstrings", type=int, nargs="+", default=[6, 8, 10, 12],
                       help="Number(s) of strings to generate layouts for")
    parser.add_argument("--spacing", type=float, default=50.0, 
                       help="Base spacing used in lattice generation")
    
    args = parser.parse_args()
    
    print("Generating lattice structure for different string counts:")
    
    for nstrings in args.nstrings:
        print(f"\n--- {nstrings} strings ---")
        x_positions, y_positions, angles = generateLatticeSpots(nstrings)
        
        print(f"Number of generated positions: {len(x_positions)}")
        print("String positions (normalized x, y, angle):")
        for i, (x, y, theta) in enumerate(zip(x_positions, y_positions, angles)):
            print(f"  String {i:2d}: ({x:6.3f}, {y:6.3f}, {theta:6.3f} rad)")
        
        print(f"Actual positions with {args.spacing}m spacing:")
        for i, (x, y) in enumerate(zip(x_positions, y_positions)):
            actual_x = x * args.spacing
            actual_y = y * args.spacing
            print(f"  String {i:2d}: ({actual_x:7.1f}, {actual_y:7.1f}) m")
            

        print("Lattice pattern (approximate):")
        grid_size = 11
        grid = [['.' for _ in range(grid_size)] for _ in range(grid_size)]
        center = grid_size // 2
            
        for i, (x, y) in enumerate(zip(x_positions, y_positions)):
            grid_x = int(round(x * 2)) + center
            grid_y = int(round(y * 2)) + center
            if 0 <= grid_x < grid_size and 0 <= grid_y < grid_size:
                grid[grid_y][grid_x] = str(i % 10)
            
        for row in grid:
            print('  ' + ' '.join(row))
