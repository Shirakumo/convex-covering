;;;; This is an implementation of the algorithm from the paper
;;;;   Convex Hull Covering of Polygonal Scenes for Accurate Collision Detection in Games
;;;; by Rong Liu et al. accessible at https://www.cs.sfu.ca/~haoz/pubs/liu_zhang_gi08.pdf
;;;;

(in-package #:org.shirakumo.fraf.convex-covering)

;;; `hull'
;;;
;;; Internal data structure for representing convex hulls during the
;;; decomposition computation. In contrast to the `convex-hull' type which is part of the API, this structure type
;;;
;;; 1. Stores vertex data using a concrete (double-float)
;;;    representation where `convex-hull' mirrors the component type
;;;    supplied by the caller
;;;
;;; 2. Stores additional information such as normals and bounding boxes

(defstruct (hull
            (:constructor %make-hull (vertices facets global-faces))
            (:copier NIL)
            (:predicate NIL))
  (vertices      (error "required") :type (manifolds:vertex-array manifolds:f64) :read-only T)
  (facets        (error "required") :type manifolds:face-array :read-only T)
  (facet-normals '())
  (bounding-box  NIL :type (or null cons)) ; (location . size/2) ; TODO maybe store 3d-space:region directly
  ;; Maybe essential
  (global-faces  #() :type manifolds:face-array :read-only T)
  ;; Debugging
  (annotations   '() :type list)
  (problem       NIL))

(defun make-hull (vertices faces vertex-position-index)
  (let ((global-faces (make-array 0 :element-type 'manifolds:u32
                                    :adjustable T
                                    :fill-pointer 0)) ; TODO remove adjust-ability before returning?
        )
    (manifolds:do-faces (a/li b/li c/li faces)
      (let* ((a    (manifolds:v vertices a/li))
             (a/gi (vertex-position a vertex-position-index))
             (b    (manifolds:v vertices b/li))
             (b/gi (vertex-position b vertex-position-index))
             (c    (manifolds:v vertices c/li))
             (c/gi (vertex-position c vertex-position-index)))
        (when (and a/gi b/gi c/gi)
          (vector-push-extend a/gi global-faces)
          (vector-push-extend b/gi global-faces)
          (vector-push-extend c/gi global-faces))))
    (%make-hull vertices faces (make-array (length global-faces)
                                           :element-type 'manifolds:u32
                                           :initial-contents global-faces))))

(defun hull-flat-p (hull &key (threshold .0001))
  (let ((vertices (hull-vertices hull)))
    (or (<= (length vertices) (* 3 3))
        (let* ((a      (manifolds:v vertices 0))
               (b      (manifolds:v vertices 1))
               (c      (manifolds:v vertices 2))
               (normal (vunit (vc (v- b a) (v- c a)))))
          (loop for i from 3 below (/ (length vertices) 3)
                for v = (manifolds:v vertices i)
                always (<= (abs (v. normal (v- v a))) threshold))))))

(defun compute-facet-normals (hull) ; TODO(jmoringe) these are just face normals for now
  (loop with vertices = (hull-vertices hull)
        with faces = (hull-facets hull)
        for index below (/ (length faces) 3)
        for normal = (vunit (manifolds:face-normal vertices faces index))
        for centroid = (manifolds:centroid vertices (subseq faces (* 3 index) (+ (* 3 index) 3)))
        ;; :when *debug-visualizations*
        ;;   :do (push (cons centroid :facet-centroid) (problem hull))
        ;;       (debug-line centroid (v* normal .3) :facet-normal hull)
        collect (cons centroid normal)))

(declaim (inline facet-normals))
(defun facet-normals (hull)
  (or (hull-facet-normals hull)
      (setf (hull-facet-normals hull) (compute-facet-normals hull))))

(defun compute-bounding-box (hull)
  (check-type hull hull)
  (let* ((vertices (hull-vertices hull))
         (faces (hull-facets hull))
         (min #1=(dvec (aref vertices 0) (aref vertices 1) (aref vertices 2)))
         (max #1#))
    (declare (type dvec3 min max))
    (manifolds:do-faces (a b c faces)
      (let ((v1 (manifolds:v vertices a))
            (v2 (manifolds:v vertices b))
            (v3 (manifolds:v vertices c)))
        (nvmin min v1 v2 v3)
        (nvmax max v1 v2 v3)))
    (center-and-size-from-min-and-max min max)))

(declaim (inline bounding-box))
(defun bounding-box (hull)
  (or (hull-bounding-box hull)
      (setf (hull-bounding-box hull)
            (multiple-value-call #'cons (compute-bounding-box hull)))))

(defun facet-bounding-box (hull facet-index)
  (check-type hull hull)
  (face-bounding-box (hull-vertices hull) (hull-facets hull) facet-index))

;;; Interface for 3d-spaces containers

(defmethod space:location ((object hull))
  (car (bounding-box object)))

(defmethod space:bsize ((object hull))
  (cdr (bounding-box object)))

;;; patch

(defstruct (patch
            (:constructor %make-patch (faces &optional surface-area hull))
            (:copier NIL)
            (:predicate NIL))
  (faces #() :type manifolds:face-array :read-only T)
  (hull NIL :type (or null hull))
  (links (make-array 0 :adjustable T :fill-pointer T) :type vector)
  (surface-area (error "required") :type double-float :read-only T)
  (compactness 0.0d0 :type double-float)
  ;; Debug information
  (link))

(defun make-patch (all-vertices a b c)
  (let* ((faces        (make-array 3 :element-type     'manifolds:u32
                                     :initial-contents (list a b c)))
         (surface-area (manifolds:face-area all-vertices faces 0)))
    (%make-patch faces surface-area)))

(defun compute-patch-convex-hull (all-vertices faces vertex-position-index)
  (let* ((vertex-count (length faces))
         ;; Quickhull doesn't like duplicate vertices so we take case
         ;; of those here.
         (vertices (make-array 3 :element-type 'manifolds:f64 :adjustable T :fill-pointer 0))
         (global-faces (make-array 3 :element-type 'manifolds:u32 :adjustable T :fill-pointer 0))
         (seen (make-hash-table :test #'equal)))
    (loop for i below vertex-count
          for j = (aref faces i)
          for a = (+ (* 3 j) 0) ; global vertex indices
          for b = (+ (* 3 j) 1)
          for c = (+ (* 3 j) 2)
          for x = (aref all-vertices a) ; TODO can we use manifold:v?
          for y = (aref all-vertices b)
          for z = (aref all-vertices c)
          for key = (list x y z)
          do (unless (gethash key seen)
               (setf (gethash key seen) T)
               #+assertions (assert (loop :for k :below (/ (length all-vertices) 3)
                                          :for v = (manifolds:v all-vertices k)
                                          :thereis (v= v (dvec x y z))))
               (vector-push-extend x vertices)
               (vector-push-extend y vertices)
               (vector-push-extend z vertices))
             (vector-push-extend j global-faces))
    (multiple-value-call #'make-hull
      (org.shirakumo.fraf.quickhull:convex-hull vertices)
      vertex-position-index)))

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
         (faces (make-array (length faces1) :element-type 'manifolds:u32 :adjustable T :fill-pointer 0))
         (seen (make-hash-table :test #'equal)))
    (flet ((add-face (a b c) ; TODO(jmoringe) check whether duplicates can actually occur
             (let ((key (list a b c)))
               (unless (gethash key seen)
                 (setf (gethash key seen) T)
                 (vector-push-extend a faces)
                 (vector-push-extend b faces)
                 (vector-push-extend c faces)))))
      (manifolds:do-faces (a b c faces1) (add-face a b c))
      (manifolds:do-faces (a b c faces2) (add-face a b c))
      (d "; Merged ~d + ~d = ~d face~:p into ~d~%"
         (length faces1) (length faces2)
         (+ (length faces1) (length faces2))
         (length faces)))
    (setf faces (make-array (length faces) :element-type 'manifolds:u32
                                           :initial-contents faces))
    (let ((all-vertices (context-vertices context)))
      (unless (= (length faces) (+ (length faces1) (length faces2))) ; TODO(jmoringe) still needed?
        (setf surface-area (manifolds:surface-area all-vertices faces)))
      (let* ((hull (compute-patch-convex-hull all-vertices faces (context-vertex-position-index context)))
             (patch (%make-patch faces surface-area hull)))
        (cond ((valid-patch-p patch context)
               (setf (patch-compactness patch) (compute-compactness all-vertices patch))
               patch)
              (T
               NIL ; hull ; TODO(jmoringe): for debugging
               ))))))

;;; `patch-link' structure
;;;
;;; A patch link references two `patch' instances that could be merged
;;; based on the mesh geometry alone. If the resulting merged patch is
;;; also valid according to the mergability criteria, that path is
;;; stored in the patch link. Finally, a patch link contains a merge
;;; cost value that is "infinite" when the merged patch is
;;; invalid. The overall algorithm performs merge operations among
;;; merge candidate path links in the order from "cheaper" candidates
;;; to "more expensive" candidates.

(defstruct (patch-link
            (:constructor %make-patch-link (a b &optional merge-cost merge-result))
            (:copier NIL)
            (:predicate NIL))
  (a NIL :type patch :read-only T)
  (b NIL :type patch :read-only T)
  (merge-cost most-positive-double-float :type double-float :read-only T)
  (merge-result NIL :type (or null patch) :read-only T))

(defun make-patch-link (context patch1 patch2)
  (assert (not (eq patch1 patch2)))
  (let ((result (make-merged-patch context patch1 patch2)))
    (if result ; (typep result 'patch)
        (let ((link (%make-patch-link patch1 patch2 (/ (patch-compactness result)) result)))
          (setf (patch-link result) link)
          link)
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

(defun merge-patches (context link)
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
                                  (make-array (length temp) :initial-contents temp :adjustable T :fill-pointer T)))
                          (unless (find-link-involving new-patch patch2)
                            (push-link new-link patch2))
                          #+no (let ((j (position link (patch-links patch2))))
                                 (setf (aref (patch-links patch2) j) new-link)))))))
      (link (patch-link-a link) (patch-link-b link))
      (link (patch-link-b link) (patch-link-a link))
      new-patch)))

(defun patch-debug-name (patch &key (format "~{[~d ~d ~d]~^ ~}~@[ …~]"))
  (let* ((faces (patch-faces patch))
         (count (length faces)))
    (format NIL format
            (coerce (subseq faces 0 (min 9 count)) 'list)
            (> count 9))))

;;; Patch link merge queue

(defun merge-priority (link)
  (let ((cost (patch-link-merge-cost link)))
    (if (= cost most-positive-double-float)
        (ash 1 31)
        (the (unsigned-byte 32) (floor cost 1/10000000)))))

(defvar *winner*)  ; the link that was selecting for merging in the current step; for visualization
(defun next-link (queue links)
  (let ((new-winner (loop for link = (damn-fast-priority-queue:dequeue queue)
                          until (or (null link)
                                    (and (gethash link links)
                                         (< (patch-link-merge-cost link) most-positive-double-float)))
                          finally (return link)))
        #+old (old-winner (loop with min = most-positive-double-float
                          with min-link = NIL
                          for link being the hash-keys of links
                          for cost = (patch-link-merge-cost link)
                          do (when (< cost min)
                               (setf min (patch-link-merge-cost link))
                               (setf min-link link))
                          finally #+no (when min-link
                                         (format *trace-output* "=> ~5,2F ~A -- ~A => ~A~%"
                                                 (patch-link-merge-cost min-link)
                                                 (patch-debug-name (patch-link-a min-link))
                                                 (patch-debug-name (patch-link-b min-link))
                                                 (patch-debug-name (patch-link-merge-result min-link))))
                                  (setf *winner* min-link)
                                  (return min-link))))
    #+old (unless (eq old-winner new-winner)
      (break "~A" (cons old-winner new-winner)))
    new-winner))

(defun link-other-patch (link this-patch)
  (cond ((eq this-patch (patch-link-a link))
         (patch-link-b link))
        ((eq this-patch (patch-link-b link))
         (patch-link-a link))
        (T
         (error "corrupt link"))))

;;; Interface

(defstruct (convex-hull
            (:constructor make-convex-hull (vertices faces))
            (:conc-name NIL)
            (:copier NIL)
            (:predicate NIL))
  (vertices (error "required") :type manifolds:vertex-array :read-only T)
  (faces    (error "required") :type manifolds:face-array :read-only T)
  (debug-info))

(defun make-convex-hull-from-hull (vertex-component-type hull)
  (let* ((hull-vertices (hull-vertices hull))
         (vertices (ecase vertex-component-type
                     (manifolds:f32
                      (map-into (make-array (length hull-vertices) :element-type 'manifolds:f32)
                                (lambda (component) (coerce component 'manifolds:f32))
                                hull-vertices))
                     (manifolds:f64
                      hull-vertices)))
         (faces (hull-facets hull)))
    (make-convex-hull vertices faces)))

(declaim (inline check-input))
(defun check-input (vertices indices)
  (check-type vertices manifolds:vertex-array)
  (check-type indices manifolds:face-array)
  (unless (zerop (mod (length indices) 3))
    (error "Total number of vertex indices in faces array is not a multiple of 3. Are all faces triangles?")))

(defun determine-vertex-component-type (vertices)
  (let ((element-type (array-element-type vertices)))
    (flet ((type= (type1 type2)
             (and (subtypep type1 type2) (subtypep type2 type1))))
      (cond ((type= element-type 'manifolds:f32)
             'manifolds:f32)
            ((type= element-type 'manifolds:f64)
             'manifolds:f64)
            (t
             ;; Should not happen since `check-input' checks the
             ;; vertex array.
             (error "Unsupported vertex component type ~s" element-type))))))

(defun coerce-input (vertices vertex-component-type)
  (ecase vertex-component-type
    (manifolds:f32
     (map-into (make-array (length vertices) :element-type 'manifolds:f64)
               (lambda (component) (coerce component 'manifolds:f64))
               vertices))
    (manifolds:f64
     vertices)))

(defun decompose (vertices indices &key) ; TODO(jmoringe) indices -> faces
  (check-input vertices indices)
  ;; FIXME: This is all really dumb and uses really bad data structures
  ;;        Could definitely be optimised a lot by someone smarter
  (let* ((vertex-component-type (determine-vertex-component-type vertices))
         (vertices (coerce-input vertices vertex-component-type))
         (context (make-context vertices indices))
         (patches (make-hash-table :test 'eq))
         (links (make-hash-table :test 'eq))

         (merge-queue (damn-fast-priority-queue:make-queue))

         (*winner* NIL))
    ;; 1. Destructure the mesh into one patch per face
    (let ((patchlist (make-array (truncate (length indices) 3)))
          (i 0))
      (manifolds:do-faces (a b c indices)
        (let ((patch (make-patch vertices a b c)))
                                        ; (setf (patch-hull patch) (compute-patch-convex-hull vertices (patch-faces patch)))
          (setf (gethash patch patches) T)
          (setf (aref patchlist i) patch)
          (incf i)))
      ;; 2. Find neighbouring patches and create the links
      (let ((*debug-visualizations* NIL))
        (let ((adjacents (manifolds:face-adjacency-list indices)))
          (dotimes (face (length adjacents))
            (loop for other in (aref adjacents face)
                  for link = (link-patches context
                                           (aref patchlist face) (aref patchlist other))
                  when link
                    do (setf (gethash link links) T)
                       (damn-fast-priority-queue:enqueue merge-queue link (merge-priority link)))))))
    ;; (visualize-step patches 0)
    ;; 3. Greedily merge patches according to merge cost
    (let ((i 1))
      (unwind-protect
           (loop for link = (progn
                              (when (zerop (mod i 100))
                                (format *trace-output* "--------Step ~:d | ~:d patch~:p ~:d link~:p~%"
                                        i (hash-table-count patches) (hash-table-count links)))
                              (next-link merge-queue links))
                 while link   ; TODO for after while is not conforming
                 for patch1 = (patch-link-a link)
                 for patch2 = (patch-link-b link)
                 for patch = (let ((*debug-visualizations* (debug-visualizations-p i)))
                               (merge-patches context link))


                 do (when (debug-visualizations-p i)
                       ;; (valid-patch-p (patch-link-merge-result link) vertices indices)
                       ;; (visualize-step patches i :highlight patch)
                       )
                     #+no (when (= i 10)
                            (loop for link in (alexandria:hash-table-keys links)
                                  for j from 0
                                  unless (typep (patch-link-merge-result link) 'patch)
                                  do (visualize-problem link i j)))

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
                    (d "  after removing ~:d patch~:p ~d link~:p~%"
                       (hash-table-count patches) (hash-table-count links))
                    ;; 2. Merge the patches
                    (let (#+no (patch (let ((*debug-visualizations* (debug-visualizations-p i)))
                                        (merge-patches vertices indices vertex-index link))))
                      ;; 3. Insert the new links
                      (assert (patch-hull patch))
                      (setf (gethash patch patches) T)
                      (loop for link across (patch-links patch)
                            do (setf (gethash link links) T)
                               (damn-fast-priority-queue:enqueue merge-queue link (merge-priority link))))
                 do                    ; consistency check
                    #+no (let ((linked-patches (make-hash-table :test #'eq))
                               (seen (make-hash-table :test #'eq)))
                           ;; validate `links' vs patch-links for all patches
                           (loop with worklist = (alexandria:hash-table-keys links)
                                 for link = (pop worklist)
                                 while link
                                 do (unless (gethash link seen)
                                      (setf (gethash link seen) T)
                                      (setf (gethash (patch-link-a link) linked-patches) T
                                            (gethash (patch-link-b link) linked-patches) T)
                                      (setf worklist (nconc worklist
                                                            (coerce (patch-links (patch-link-a link)) 'list)
                                                            (coerce (patch-links (patch-link-b link)) 'list)))))
                           (unless (alexandria:set-equal (alexandria:hash-table-keys patches)
                                                         (alexandria:hash-table-keys linked-patches))
                             (format *error-output* "~d patches ~d linked patches~%"
                                     (hash-table-count patches)
                                     (hash-table-count linked-patches)))
                           (unless (= (hash-table-count patches) 1)
                             (assert (alexandria:set-equal (alexandria:hash-table-keys patches)
                                                           (alexandria:hash-table-keys linked-patches)))))
                    (incf i)
                 finally (format *trace-output* "--------Result | ~:d patch~:p ~:d link~:p~%"
                                  (hash-table-count patches) (hash-table-count links)))
        ;; (visualize-step (alexandria:hash-table-keys patches) i :final T)
        ))
    ;; 4. Return the patches' convex hulls
    (loop with hulls = (make-array (hash-table-count patches))
          for i from 0
          for patch being the hash-keys of patches
          for hull = (patch-hull patch)
          do (setf (aref hulls i) (when hull ; TODO(jmoringe): can we avoid storing those patches in the first place?
                                    (let ((result (make-convex-hull-from-hull vertex-component-type hull)))
                                      (setf (debug-info result) (list :hull hull :patch patch))
                                      result)))
          finally (return (values (remove nil hulls) context)))))
