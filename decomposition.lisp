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
            (:constructor %make-patch (faces &optional hull compactness links))
            (:copier NIL)
            (:predicate NIL))
  (faces #() :type (simple-array (unsigned-byte 32) (*)))
  (hull NIL :type (or null convex-hull))
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

(defun compute-compactness (vertices patch)
  (/ (sqrt (manifolds:surface-area vertices (patch-faces patch)))
     (manifolds:boundary-length vertices (patch-faces patch))))

(defun make-merged-patch (vertices a b)
  (let ((faces (make-array (+ (length (patch-faces a)) (length (patch-faces a))) :element-type '(unsigned-byte 32))))
    (replace faces (patch-faces a))
    (replace faces (patch-faces b) :start1 (length (patch-faces a)))
    (multiple-value-bind (hull-verts hull-faces) (org.shirakumo.fraf.quickhull:convex-hull vertices faces)
      (let ((patch (%make-patch faces (make-convex-hull hull-verts hull-faces))))
        (when (valid-patch-p patch)
          (setf (patch-compactness patch) (compute-compactness vertices patch))
          patch)))))

(defun make-patch (a b c)
  (%make-patch (make-array 3 :element-type '(unsigned-byte 32) :initial-contents (list a b c))))

(defstruct (patch-link
            (:constructor %make-patch-link (a b &optional merge-cost merge-result))
            (:copier NIL)
            (:predicate NIL))
  (a NIL :type patch)
  (b NIL :type patch)
  (merge-cost most-positive-double-float :type double-float)
  (merge-result NIL :type (or null patch)))

(defun make-patch-link (vertices a b)
  (let ((result (make-merged-patch vertices a b)))
    (if result
        (%make-patch-link a b (/ (patch-compactness result)) result)
        (%make-patch-link a b))))

(defun link-patches (vertices a b)
  (unless (loop for link across (patch-links a)
                thereis (or (eql b (patch-link-a link)) (eql b (patch-link-b link))))
    (let ((link (make-patch-link vertices a b)))
      (vector-push-extend link (patch-links a))
      (vector-push-extend link (patch-links b)))))

(defun merge-patches (vertices link)
  (let ((patch (patch-link-merge-result link)))
    ;; Now that we actually merge this in, compute new links and update the existing neighbour's links.
    (flet ((link (cur other)
             (loop for i from 0 below (length (patch-links cur))
                   for link = (aref (patch-links cur) i)
                   do (unless (or (eql other (patch-link-a link)) (eql other (patch-link-b link)))
                        (let ((new-link (make-patch-link vertices patch (if (eql cur (patch-link-a link)) (patch-link-b link) (patch-link-a link)))))
                          (vector-push-extend new-link (patch-links patch))
                          (setf (aref (patch-links cur) i) new-link))))))
      (link (patch-link-a link) (patch-link-b link))
      (link (patch-link-b link) (patch-link-a link))
      patch)))

(defun next-link (links)
  (loop with min = most-positive-double-float
        with min-link = NIL
        for link being the hash-keys of links
        do (when (<= (patch-link-merge-cost link) min)
             (setf min (patch-link-merge-cost link))
             (setf min-link link))
        finally (return min-link)))

(defun decompose (vertices indices)
  ;; FIXME: This is all really dumb and uses really bad data structures
  ;;        Could definitely be optimised a lot by someone smarter
  (let ((patchlist (make-array (truncate (length indices) 3)))
        (patches (make-hash-table :test 'eq))
        (links (make-hash-table :test 'eq))
        (i 0))
    ;; 1. Destructure the mesh into one patch per face
    (manifolds:do-faces (a b c indices)
      (let ((patch (make-patch a b c)))
        (setf (gethash patch patches) T)
        (setf (aref patchlist i) patch)
        (incf i)))
    ;; 2. Find neighbouring patches and create the links
    (let ((adjacents (manifolds:face-adjacency-list indices)))
      (dotimes (face (length adjacents))
        (loop for other in (aref adjacents face)
              do (link-patches vertices (aref patchlist face) (aref patchlist other)))))
    ;; 3. Greedily merge patches according to merge cost
    (loop for link = (next-link links)
          while (patch-link-merge-result link)
          do ;; 1. Remove old links
          (remhash link links)
          (remhash (patch-link-a link) patches)
          (remhash (patch-link-b link) patches)
          (loop for other across (patch-links (patch-link-b link))
                do (remhash other links))
          (loop for other across (patch-links (patch-link-a link))
                do (remhash other links))
          ;; 2. Merge the patches
          (let ((patch (merge-patches vertices link)))
            ;; 3. Insert the new links
            (setf (gethash patch patches) T)
            (loop for link across (patch-links patch)
                  do (setf (gethash link links) T))))
    ;; 4. Return the patches' convex hulls
    (let ((hulls (make-array (hash-table-count patches))))
      (loop for patch being the hash-keys of patches
            for i from 0
            do (setf (aref hulls i) (patch-hull patch)))
      hulls)))
