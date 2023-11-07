;;;; This is an implementation of the algorithm from the paper
;;;;   Convex Hull Covering of Polygonal Scenes for Accurate Collision Detection in Games
;;;; by Rong Liu et al. accessible at https://www.cs.sfu.ca/~haoz/pubs/liu_zhang_gi08.pdf
;;;;

(in-package #:org.shirakumo.fraf.convex-covering)

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
