;;;; This is an implementation of the algorithm from the paper
;;;;   Convex Hull Covering of Polygonal Scenes for Accurate Collision Detection in Games
;;;; by Rong Liu et al. accessible at https://www.cs.sfu.ca/~haoz/pubs/liu_zhang_gi08.pdf
;;;;

(in-package #:org.shirakumo.fraf.convex-covering)

;;; Vertex index
;;;
;;; Lookup index in global vertex array given the vertex coordinates.

(defun make-vertex-index ()
  (make-hash-table :test #'equalp))

(defun vertex-position (vertex index)
  (check-type vertex dvec3)
  (the (or null manifolds:u32) (gethash vertex index)))

(defun (setf vertex-position) (new-value vertex index)
  (check-type vertex dvec3)
  (setf (gethash vertex index) new-value))

(defun index-vertices (vertices)
  (check-type vertices manifolds:vertex-array)
  (loop :with index = (make-vertex-index)
        :for i :below (/ (length vertices) 3)
        :for vertex = (manifolds:v vertices i)
        :do (setf (vertex-position vertex index) i)
        :finally (return index)))

;;; Face index

(defstruct (face-info
            (:constructor %make-face-info (index center size/2)))
  (index  (error "required") :type manifolds:u32)
  (center (error "required") :type vec3)
  (size/2 (error "required") :type vec3))

(defun make-face-info (face-index v1 v2 v3)
  (declare (type dvec3 v1 v2 v3))
  (let* ((min    (vmin v1 v2 v3))
         (max    (vmax v1 v2 v3))
         (center (v/ (v+ min max) 2))
         (size/2 (v- center min)))
    (%make-face-info face-index
                     (vec (vx center) (vy center) (vz center))
                     (vec (vx size/2) (vy size/2) (vz size/2)))))

(defmethod space:location ((object face-info))
  (face-info-center object))

(defmethod space:bsize ((object face-info))
  (face-info-size/2 object))

(defun index-faces (vertices faces)
  (let ((index (org.shirakumo.fraf.trial.space.grid3:make-grid .3 :bsize (vec3 2 2 2))))
    (loop :for i :below (/ (length faces) 3)
          :for v1 = (manifolds:v vertices (aref faces (+ (* 3 i) 0)))
          :for v2 = (manifolds:v vertices (aref faces (+ (* 3 i) 1)))
          :for v3 = (manifolds:v vertices (aref faces (+ (* 3 i) 2)))
          :for info = (make-face-info i v1 v2 v3)
          :do (space:enter info index))
    ; (space:reoptimize index)
    (space:do-overlapping (face-info index (space:region .5 .0 .0 .9 .9 .9))
      (format *trace-output* "~A~%" face-info))
    (break "~A" index)
    index))

;;; Context

(defstruct (context
            (:constructor make-context (vertices faces
                                        &aux (vertex-index (index-vertices vertices))
                                             (face-index   (index-faces vertices faces))))
            (:copier nil)
            (:predicate nil))
  ;; Mesh
  (vertices     (error "required"))
  (faces        (error "required"))
  ;; Index structures
  (vertex-index (error "required"))
  (face-index   (error "required")))

;;;

(defstruct (convex-hull
            (:constructor %make-convex-hull (vertices faces global-faces))
            (:conc-name NIL)
            (:copier NIL)
            (:predicate NIL))
  (vertices #() :type (manifolds:vertex-array double-float))
  (faces    #() :type manifolds:face-array)
  (face-normals '())                    ; TODO facet normals
  ;; Maybe essential
  (global-faces #() :type manifolds:face-array)
  ;; Debugging
  (annotations '() :type list)
  (problem     nil))

(defun make-convex-hull (vertices faces vertex-index)
  (let ((global-faces (make-array 0 :element-type 'manifolds:u32
                                    :adjustable t
                                    :fill-pointer 0)) ; TODO remove adjust-ability before returning?
        )
    (manifolds:do-faces (a/li b/li c/li faces)
      (let* ((a    (manifolds:v vertices a/li))
             (a/gi (vertex-position a vertex-index))
             (b    (manifolds:v vertices b/li))
             (b/gi (vertex-position b vertex-index))
             (c    (manifolds:v vertices c/li))
             (c/gi (vertex-position c vertex-index)))
        (when (and a/gi b/gi c/gi)
          (vector-push-extend a/gi global-faces)
          (vector-push-extend b/gi global-faces)
          (vector-push-extend c/gi global-faces))))
    (%make-convex-hull vertices faces (make-array (length global-faces)
                                                  :element-type 'manifolds:u32
                                                  :initial-contents global-faces))))

(defstruct (patch
            (:constructor %make-patch (faces &optional surface-area hull))
            (:copier NIL)
            (:predicate NIL))
  (faces #() :type manifolds:face-array)
  (hull NIL :type (or null convex-hull))
  (links (make-array 0 :adjustable T :fill-pointer T) :type vector)
  (surface-area (error "required") :type double-float)
  (compactness 0.0d0 :type double-float))

(defun make-patch (all-vertices a b c)
  (let* ((faces        (make-array 3 :element-type     'manifolds:u32
                                     :initial-contents (list a b c)))
         (surface-area (manifolds:face-area all-vertices faces 0)))
    (%make-patch faces surface-area)))

(defun compute-patch-convex-hull (all-vertices faces vertex-index)
  (let* ((vertex-count (length faces))
         (vertices ; (make-array (* 3 vertex-count) :element-type 'double-float)
           (make-array 3 :element-type 'double-float :adjustable t :fill-pointer 0))
         (global-faces
           (make-array 3 :element-type '(unsigned-byte 32) :adjustable t :fill-pointer 0)))
    (let ((seen (make-hash-table :test #'equal)))
      (loop :for i :below vertex-count
            :for j = (aref faces i)
            :for a = (+ (* 3 j) 0) ; global vertex indices
            :for b = (+ (* 3 j) 1)
            :for c = (+ (* 3 j) 2)
            :for x = (aref all-vertices a) ; TODO can we use manifold:v?
            :for y = (aref all-vertices b)
            :for z = (aref all-vertices c)
            :for key = (list x y z)
            :do (unless (gethash key seen) ; in case quickhull does not permit duplicate vertices
                  (setf (gethash key seen) t)
                  (assert (loop :for k :below (/ (length all-vertices) 3)
                                :for v = (manifolds:v all-vertices k)
                                :thereis (v= v (dvec x y z))))
                  (vector-push-extend x vertices)
                  (vector-push-extend y vertices)
                  (vector-push-extend z vertices))
                (vector-push-extend j global-faces)))
    (multiple-value-call #'make-convex-hull
      (org.shirakumo.fraf.quickhull:convex-hull vertices)
      vertex-index)))

(defun push-link (new-link patch)
  (assert (or (eq patch (patch-link-a new-link))
              (eq patch (patch-link-b new-link))))
  (assert (not (find-if (lambda (link)
                          (or (and (eq (patch-link-a link) (patch-link-a new-link))
                                   (eq (patch-link-b link) (patch-link-b new-link)))
                              (and (eq (patch-link-b link) (patch-link-a new-link))
                                   (eq (patch-link-a link) (patch-link-b new-link)))))
                        (patch-links patch))))
  (vector-push-extend new-link (patch-links patch)))

(defun find-link-involving (other-patch patch)
  (find-if (lambda (link)
             (link-involves-p other-patch link))
           (patch-links patch)))

;;;

(defun make-merged-patch (context a b)
  (assert (not (eq a b)))
  (let* ((surface-area (+ (patch-surface-area a) (patch-surface-area b))) ; TODO only when valid
         (faces1 (patch-faces a))
         (faces2 (patch-faces b))
         (faces (make-array (length faces1) :element-type 'manifolds:u32 :adjustable t :fill-pointer 0))
         (seen (make-hash-table :test #'equal)))
    (flet ((add-face (a b c) ; TODO check whether duplicates can actually occur
             (let ((key (list a b c)))
               (unless (gethash key seen)
                 (setf (gethash key seen) t)
                 (vector-push-extend a faces)
                 (vector-push-extend b faces)
                 (vector-push-extend c faces)))))
      (manifolds:do-faces (a b c faces1) (add-face a b c))
      (manifolds:do-faces (a b c faces2) (add-face a b c))
      (d "; Merged ~D + ~D = ~D face~:P into ~D~%"
         (length faces1) (length faces2)
         (+ (length faces1) (length faces2))
         (length faces)))
    (setf faces (make-array (length faces) :element-type 'manifolds:u32
                                           :initial-contents faces))
    (let ((all-vertices (context-vertices context)))
      (unless (= (length faces) (+ (length faces1) (length faces2))) ; TODO still needed?
        (setf surface-area (manifolds:surface-area all-vertices faces)))
      (let* ((hull (compute-patch-convex-hull all-vertices faces (context-vertex-index context)))
             (patch (%make-patch faces surface-area hull)))
        (cond ((valid-patch-p patch context)
               (setf (patch-compactness patch) (compute-compactness all-vertices patch))
               patch)
              (t                        ; TODO for debugging
               hull))))))

(defstruct (patch-link
            (:constructor %make-patch-link (a b &optional merge-cost merge-result))
            (:copier NIL)
            (:predicate NIL))
  (a NIL :type patch)
  (b NIL :type patch)
  (merge-cost most-positive-double-float :type double-float)
  (merge-result NIL ; :type (or null patch)
   ))

(defun make-patch-link (context patch1 patch2)
  (assert (not (eq patch1 patch2)))
  (let ((result (make-merged-patch context patch1 patch2)))
    (if (typep result 'patch)
        (%make-patch-link patch1 patch2 (/ (patch-compactness result)) result)
        (%make-patch-link patch1 patch2 most-positive-double-float result))))

(defun link-patches (context a b) ; called in preparation phase for adjacent faces
  (unless (loop for link across (patch-links a) ; TODO make a predicate
                thereis (or (eql b (patch-link-a link))
                            (eql b (patch-link-b link))))
    (let ((link (make-patch-link context a b)))
      (push-link link a)
      (push-link link b)
      link)))

(defun link-involves-p (patch link)
  (or (eq (patch-link-a link) patch)
      (eq (patch-link-b link) patch)))

(defun merge-patches (all-vertices all-faces context link)
  (let ((new-patch (patch-link-merge-result link)))
    ;; Now that we actually merge this in, compute new links and update the existing neighbour's links.
    (flet ((link (cur other)
             (loop with links = (patch-links cur)
                   for i from 0 below (length links)
                   for link = (aref links i)
                   do (unless (link-involves-p other link)
                        (let* ((patch2 (link-other-patch link cur))
                               (new-link (make-patch-link context new-patch patch2)))
                          (unless (find-link-involving patch2 new-patch)
                            (push-link new-link new-patch))
                          (let ((temp (remove-if (lambda (old-link)
                                                   ; (assert (not (link-involves-p new-patch old-link)))
                                                   (when (eq old-link link)
                                                     (assert (link-involves-p cur old-link)))
                                                   (link-involves-p cur old-link)
                                                   ;; (link-involves-p other old-link)
                                                   )
                                                 (patch-links patch2))))
                            (setf (patch-links patch2)
                                  (make-array (length temp) :initial-contents temp :adjustable t :fill-pointer t)))
                          (unless (find-link-involving new-patch patch2)
                            (push-link new-link patch2))
                          #+no (let ((j (position link (patch-links patch2))))
                                 (setf (aref (patch-links patch2) j) new-link)))))))
      (link (patch-link-a link) (patch-link-b link))
      (link (patch-link-b link) (patch-link-a link))
      new-patch)))

(defun patch-debug-name (patch &key (format "~{[~D ~D ~D]~^ ~}~@[ …~]"))
  (let* ((faces (patch-faces patch))
         (count (length faces)))
    (format nil format
            (coerce (subseq faces 0 (min 9 count)) 'list)
            (> count 9))))

(defun next-link (links)
  (loop with min = most-positive-double-float
        with min-link = NIL
        for link being the hash-keys of links
        for cost = (patch-link-merge-cost link)
        do (when (< cost min)
             (setf min (patch-link-merge-cost link))
             (setf min-link link))
        finally (when min-link
                  (format *trace-output* "=> ~5,2F ~A -- ~A => ~A~%"
                          (patch-link-merge-cost min-link)
                          (patch-debug-name (patch-link-a min-link))
                          (patch-debug-name (patch-link-b min-link))
                          (patch-debug-name (patch-link-merge-result min-link))))
                (setf *winner* min-link)
                (return min-link)))

(defun link-other-patch (link this-patch)
  (cond ((eq this-patch (patch-link-a link))
         (patch-link-b link))
        ((eq this-patch (patch-link-b link))
         (patch-link-a link))
        (t
         (error "corrupt link"))))

(defun verify-input (vertices indices)
  (unless (zerop (mod (length indices) 3))
    (error "Total number of vertex indices in faces array is not a multiple of 3. Are all faces triangles?")))

(defun decompose (vertices indices &key) ; TODO indices -> faces
  (verify-input vertices indices)
  ;; FIXME: This is all really dumb and uses really bad data structures
  ;;        Could definitely be optimised a lot by someone smarter
  (let* ((context (make-context vertices indices))
         (patchlist (make-array (truncate (length indices) 3)))
         (patches (make-hash-table :test 'eq))
         (links (make-hash-table :test 'eq))
         (i 0)
         (*winner* nil))
    ;; 1. Destructure the mesh into one patch per face
    (manifolds:do-faces (a b c indices)
      (let ((patch (make-patch vertices a b c)))
                                        ; (setf (patch-hull patch) (compute-patch-convex-hull vertices (patch-faces patch)))
        (setf (gethash patch patches) T)
        (setf (aref patchlist i) patch)
        (incf i)))
    ;; 2. Find neighbouring patches and create the links
    (let ((adjacents (manifolds:face-adjacency-list indices)))
      (dotimes (face (length adjacents))
        (loop for other in (aref adjacents face)
              for link = (link-patches context
                                       (aref patchlist face) (aref patchlist other))
              when link
              do (setf (gethash link links) T))))
    (visualize-step patches 0)
    ;; 3. Greedily merge patches according to merge cost
    (let ((i 1))
      (unwind-protect
           (loop for link = (progn
                              (format *trace-output* "--------Step ~D | ~:D patch~:P ~D link~:P~%"
                                      i (hash-table-count patches) (hash-table-count links))
                              (next-link links))
                 while link   ; TODO for after while is not conforming
                 for patch1 = (patch-link-a link)
                 for patch2 = (patch-link-b link)
                 for patch = (let ((*debug-visualizations* (debug-visualizations-p i)))
                               (merge-patches vertices indices context link))


                 :do (when (debug-visualizations-p i)
                       ;; (valid-patch-p (patch-link-merge-result link) vertices indices)
                       (visualize-step patches i :highlight patch))
                     #+no (when (= i 10)
                            (loop :for link :in (alexandria:hash-table-keys links)
                                  :for j :from 0
                                  :unless (typep (patch-link-merge-result link) 'patch)
                                  :do (visualize-problem link i j)))

                 do ;; 1. Remove old patches and links
                    (assert (not (eq patch1 patch2)))
                    (assert (remhash patch1 patches))
                    (assert (remhash patch2 patches))
                    (remhash link links)
                    (flet ((remove-other-links (patch)
                             (loop for other-link across (patch-links patch)
                                        ; for other-patch = (link-other-patch other-link patch)
                                   do (remhash other-link links)
                                        ; (setf (patch-links other-patch) (remove other-link (patch-links other-patch)))
                                   )))
                      (remove-other-links patch1)
                      (remove-other-links patch2))
                    (d "  after removing ~:D patch~:P ~D link~:P~%"
                       (hash-table-count patches) (hash-table-count links))
                    ;; 2. Merge the patches
                    (let (#+no (patch (let ((*debug-visualizations* (debug-visualizations-p i)))
                                        (merge-patches vertices indices vertex-index link))))
                      ;; 3. Insert the new links
                      (assert (patch-hull patch))
                      (setf (gethash patch patches) T)
                      (loop for link across (patch-links patch)
                            do (setf (gethash link links) T)))
                 :do                    ; consistency check
                    #+no (let ((linked-patches (make-hash-table :test #'eq))
                               (seen (make-hash-table :test #'eq)))
                           ;; validate `links' vs patch-links for all patches
                           (loop :with worklist = (alexandria:hash-table-keys links)
                                 :for link = (pop worklist)
                                 :while link
                                 :do (unless (gethash link seen)
                                       (setf (gethash link seen) t)
                                       (setf (gethash (patch-link-a link) linked-patches) t
                                             (gethash (patch-link-b link) linked-patches) t)
                                       (setf worklist (nconc worklist
                                                             (coerce (patch-links (patch-link-a link)) 'list)
                                                             (coerce (patch-links (patch-link-b link)) 'list)))))
                           (unless (alexandria:set-equal (alexandria:hash-table-keys patches)
                                                         (alexandria:hash-table-keys linked-patches))
                             (format *error-output* "~D patches ~D linked patches~%"
                                     (hash-table-count patches)
                                     (hash-table-count linked-patches)))
                           (unless (= (hash-table-count patches) 1)
                             (assert (alexandria:set-equal (alexandria:hash-table-keys patches)
                                                           (alexandria:hash-table-keys linked-patches)))))
                    (incf i)
                 :finally (format *trace-output* "--------Result | ~:D patch~:P ~D link~:P~%"
                                  (hash-table-count patches) (hash-table-count links)))
        (visualize-step (alexandria:hash-table-keys patches) i :final t)))
    ;; 4. Return the patches' convex hulls
    (let ((hulls (make-array (hash-table-count patches))))
      (loop for patch being the hash-keys of patches
            for i from 0
            do (setf (aref hulls i) (patch-hull patch)))
      (remove nil hulls)                ; TODO temp hack
      )))
