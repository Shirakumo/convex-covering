;;;; This is an implementation of the algorithm from the paper
;;;;   Convex Hull Covering of Polygonal Scenes for Accurate Collision Detection in Games
;;;; by Rong Liu et al. accessible at https://www.cs.sfu.ca/~haoz/pubs/liu_zhang_gi08.pdf
;;;;

(in-package #:org.shirakumo.fraf.convex-covering)

(defstruct (convex-hull
            (:constructor make-convex-hull (vertices faces))
            (:copier NIL)
            (:predicate NIL))
  (vertices  #() :type (simple-array double-float (*)))
  (faces #() :type (simple-array (unsigned-byte 32) (*))))

(defun find-vertex-in-hull (hull)
  )

(defun find-edge-in-hull (hull)
  )

(defun find-boundary-constraint-bar (hull)
  )

(defun find-negative-side-touching-triangle (hull)
  )

(defstruct (patch
            (:constructor %make-patch (faces hull &optional compactness))
            (:copier NIL)
            (:predicate NIL))
  (faces #() :type (simple-array (unsigned-byte 32) (*)))
  (hull NIL :type convex-hull)
  (links (make-array 0 :adjustable T :fill-pointer T) :type vector)
  (compactness 0.0d0 :type double-float))

(defun normals-matching-p (patch hull)
  )

(defun valid-patch-p (patch)
  (let ((hull (patch-hull patch)))
    (and (not (find-vertex-in-hull hull))
         (not (find-edge-in-hull hull))
         (normals-matching-p patch hull)
         (not (find-boundary-constraint-bar hull))
         (not (find-negative-side-touching-triangle hull)))))

(defun compute-compactness (patch)
  (/ (sqrt (surface-area ... (patch-faces patch)))
     (boundary-length ... (patch-faces patch))))

(defun make-merged-patch (a b)
  (let ((faces (make-array (+ (length (patch-faces a)) (length (patch-faces a))) :element-type '(unsigned-byte 32))))
    (replace faces (patch-faces a))
    (replace faces (patch-faces b) :start1 (length (patch-faces a)))
    (multiple-value-bind (hull-verts hull-faces) (org.shirakumo.fraf.quickhull:convex-hull verts faces)
      (let ((patch (%make-patch faces (%make-convex-hull hull-verts hull-faces))))
        (when (valid-patch-p patch)
          (setf (patch-compactness patch) (compute-compactness patch))
          patch)))))

(defstruct (patch-link
            (:constructor %make-patch-link (a b &optional merge-cost merge-result))
            (:copier NIL)
            (:predicate NIL))
  (a NIL :type patch)
  (b NIL :type patch)
  (merge-cost most-positive-double-float :type double-float)
  (merge-result NIL :type (or null patch)))

(defun make-patch-link (a b)
  (let ((result (make-merged-patch a b)))
    (if result
        (%make-patch-link a b (/ (patch-compactness result)) result)
        (%make-patch-link a b))))

(defun merge-patches (link)
  (let ((patch (patch-link-merge-result link)))
    ;; Now that we actually merge this in, compute new links and update the existing neighbour's links.
    (flet ((link (cur other)
             (loop for i from 0 below (length (patch-links cur))
                   for link = (aref (patch-links cur) i)
                   do (unless (or (eql other (patch-link-a link)) (eql other (patch-link-b link)))
                        (let ((new-link (make-patch-link patch (if (eql cur (patch-link-a link)) (patch-link-b link) (patch-link-a link)))))
                          (vector-push-extend new-link (patch-links patch))
                          (setf (aref (patch-links cur) i) new-link))))))
      (link (patch-link-a link) (patch-link-b link))
      (link (patch-link-b link) (patch-link-a link))
      patch)))

(defun decompose (vertices indices)
  ;; 1. Destructure the mesh into one patch per face
  ;; 2. Find neighbouring patches and create the links
  ;; 3. Greedily merge patches according to merge cost
  ;; 4. Return the patches' convex hulls
  )
