import rhinoscriptsyntax as rs
import scriptcontext as sc
import Rhino

def get_face_center(face):
    bbox = face.GetBoundingBox(True)
    if not bbox:
        return None
    return bbox.Center

def get_lower_profile_and_export():
    # Select one extrusion or polysurface
    obj = rs.GetObject("Select extrusion or solid", rs.filter.polysurface | rs.filter.extrusion)
    if not obj:
        return

    brep = rs.coercebrep(obj)
    if not brep:
        print("Not a Brep or extrusion.")
        return

    # Find all planar faces and their Z center
    planar_faces = []
    for face in brep.Faces:
        if face.IsPlanar():
            center = get_face_center(face)
            if center:
                planar_faces.append((center.Z, face))

    if not planar_faces:
        print("No planar faces found.")
        return

    # Sort to find lowest face
    planar_faces.sort(key=lambda x: x[0])
    lowest_face = planar_faces[0][1]

    # Get outer loop of the face
    loop = lowest_face.OuterLoop
    curve = loop.To3dCurve()
    if not curve:
        print("Could not extract curve from outer loop.")
        return

    # Add curve to Rhino (optional visualization)
    sc.doc.Objects.AddCurve(curve)
    sc.doc.Views.Redraw()
    print("Lower profile extracted.")

    # Ask user for file to export
    filepath = rs.SaveFileName("Save GSE File", "GSE Files (*.gse)|*.gse||")
    if not filepath:
        return

    export_profile_to_gse_approx(curve, filepath)
    print("Exported to {}".format(filepath))

def export_profile_to_gse_approx(curve, filepath, segment_count=50):
    # Divide curve into points (approximate with segments)
    points = curve.DivideByCount(segment_count, True)
    if not points:
        print("Failed to divide curve.")
        return

    pt_list_3d = [curve.PointAt(t) for t in points]

    # Get curve's underlying plane (the curve is planar, so this works)
    rc, plane = curve.TryGetPlane()
    if not rc:
        print("Curve is not planar or can't get plane.")
        return

    # Convert 3D points to 2D plane coordinates (U,V)
    pt_list_2d = [plane.ClosestPoint(pt) for pt in pt_list_3d]  # Returns Point3d but U,V in X,Y components

    with open(filepath, "w") as f:
        f.write("$debut du fichier\n")
        f.write("$version\n1\n")
        f.write("$ulong\nunit\n")  # Optional: you can put your units here
        f.write("$points\n")
        f.write("{}\n".format(len(pt_list_2d)))

        # Export U and V (local plane X and Y coords)
        for i, pt in enumerate(pt_list_2d):
            # pt.X = U, pt.Y = V
            f.write("{} {} {}\n".format(i + 1, pt.X, pt.Y))

        f.write("$courbes\n")
        for i in range(len(pt_list_2d) - 1):
            f.write("segment {} {}\n".format(i + 1, i + 2))
        if rs.IsCurveClosed(curve):
            f.write("segment {} {}\n".format(len(pt_list_2d), 1))

        f.write("$fin du fichier\n")

    print("Exported profile with local plane coordinates to {}".format(filepath))


# Run the full process
get_lower_profile_and_export()
