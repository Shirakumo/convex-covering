;;;; This is an implementation of the algorithm from the paper
;;;;   Convex Hull Covering of Polygonal Scenes for Accurate Collision Detection in Games
;;;; by Rong Liu et al. accessible at https://www.cs.sfu.ca/~haoz/pubs/liu_zhang_gi08.pdf
;;;;

(in-package #:org.shirakumo.fraf.convex-covering)

;;; Patch link merge queue

(defun merge-priority (link)
  (let ((cost (patch-link-merge-cost link)))
    (unless (= cost most-positive-double-float)
      ;; The multiplication factor 10000000 has been determined
      ;; empirically. The MIN ensures that the result will actually
      ;; fit into 32 bits. Taking the MIN will of course squash the
      ;; differences between certain, high costs but those should end
      ;; up at the end of the queue and thus not matter much .
      (the (unsigned-byte 32)
           (values (floor (min cost
                               #.(coerce (/ #.(ash 1 30) #1=10000000)
                                         'double-float))
                          #.(/ 1 #1#)))))))

(defun maybe-enqueue-link (link queue)
  (let ((priority (merge-priority link)))
    (when priority
      (damn-fast-priority-queue:enqueue queue link priority))))

(defun next-link (queue links)
  (loop for link = (damn-fast-priority-queue:dequeue queue)
        until (or (null link)
                  (and (gethash link links)
                       (< (patch-link-merge-cost link) most-positive-double-float))) ; TODO(jmoringe): should not happen
        finally (return link)))

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
  (check-type indices vector)
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

(defun coerce-input (vertices faces vertex-component-type normalization-threshold)
  (let ((vertices (ecase vertex-component-type
                    (manifolds:f32
                     (map-into (make-array (length vertices) :element-type 'manifolds:f64)
                               (lambda (component) (coerce component 'manifolds:f64))
                               vertices))
                    (manifolds:f64
                     vertices)))
        (faces (coerce faces 'manifolds:face-array)))
    (manifolds:normalize vertices faces :threshold normalization-threshold)))

(defun decompose (vertices faces &key (merge-cost              #'/compactness)
                                      (tolerance               1d-4)
                                      (normalization-threshold (/ tolerance 10d0))
                                      (edge-tolerance          tolerance)
                                      (normals-tolerance       1d-4 #+TODO tolerance)
                                      (thread-pool             T))
  (check-input vertices faces)
  (let ((vertex-component-type (determine-vertex-component-type vertices)))
    (multiple-value-bind (vertices faces)
        (coerce-input vertices faces vertex-component-type normalization-threshold)
      (with-thread-pool (thread-pool)
        (let* ((context (make-context vertices faces (coerce-to-cost-function merge-cost)
                                      :edge-tolerance    (coerce edge-tolerance 'double-float)
                                      :normals-tolerance (coerce normals-tolerance 'double-float)))
               (patches (make-hash-table :test 'eq))
               (links (make-hash-table :test 'eq))
               (merge-queue (damn-fast-priority-queue:make-queue)))
          ;; 1. Destructure the mesh into one patch per face
          (let ((patchlist (make-array (truncate (length faces) 3)))
                (i 0))
            (manifolds:do-faces (a b c faces)
              (let ((patch (make-initial-patch vertices a b c)))
                (setf (gethash patch patches) T)
                (setf (aref patchlist i) patch)
                (incf i)))
            ;; 2. Find neighbouring patches and create the links
            (let ((*debug-visualizations* NIL))
              (let ((adjacents (manifolds:face-adjacency-list faces)))
                (dotimes (face (length adjacents))
                  (loop for other in (aref adjacents face)
                        for link = (link-patches context
                                                 (aref patchlist face) (aref patchlist other))
                        when link
                          do (setf (gethash link links) T) ; TODO(jmoringe): don't store if invalid?
                             (let ((priority (merge-priority link)))
                               (when priority
                                 (damn-fast-priority-queue:enqueue merge-queue link priority))))))))
          ;; 3. Greedily merge patches according to merge cost
          (flet ((report (step)
                   (format *trace-output* "-------- ~:[Result~;Step ~:*~:d~] | ~:d patch~:p ~:d link~:p~%"
                           step (hash-table-count patches) (hash-table-count links))))
            (loop with i = 1
                  for link = (progn
                               (when (zerop (mod i 100)) (report i))
                               (next-link merge-queue links))
                  while link  ; TODO for after while is not conforming
                  for patch1 = (patch-link-a link)
                  for patch2 = (patch-link-b link)
                  ;; a) Merge the patches
                  for patch = (let ((*debug-visualizations* (debug-visualizations-p i)))
                                (merge-patches context link))
                  ;; b) Remove old patches and links
                  do (assert (not (eq patch1 patch2)))
                     (assert (remhash patch1 patches))
                     (assert (remhash patch2 patches))
                     (remhash link links)
                     (flet ((remove-other-links (patch)
                              (loop for other-link across (patch-links patch)
                                    do (remhash other-link links))))
                       (remove-other-links patch1)
                       (remove-other-links patch2))
                     (d "  after removing ~:d patch~:p ~d link~:p~%"
                        (hash-table-count patches) (hash-table-count links))
                     ;; c) Insert the new links
                     (assert (patch-hull patch))
                     (setf (gethash patch patches) T)
                     (loop for link across (patch-links patch)
                           do (setf (gethash link links) T)
                              (let ((priority (merge-priority link)))
                                (when priority
                                  (damn-fast-priority-queue:enqueue merge-queue link priority))))
                  do ; consistency check
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
                  finally (report nil)))
          ;; 4. Return the patches' convex hulls
          (loop with hulls = (make-array (hash-table-count patches))
                for i from 0
                for patch being the hash-keys of patches
                for hull = (patch-hull patch)
                ;; A patch may not have a hull if the patch never got
                ;; merged. Compute those missing hulls now.
                do (setf (aref hulls i)
                         (handler-case
                             (let ((hull (or hull
                                             (compute-patch-convex-hull
                                              vertices
                                              (patch-faces patch)
                                              (context-vertex-position-index context)))))
                               (let ((result (make-convex-hull-from-hull vertex-component-type hull)))
                                 (setf (debug-info result) (list :hull hull :patch patch))
                                 result))
                           (org.shirakumo.fraf.quickhull:points-colinear-error (condition)
                             (warn "Bad hull: ~A" condition)
                             nil)))
                finally (return (values (remove nil hulls) context))))))))
