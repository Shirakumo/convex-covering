(in-package #:org.shirakumo.fraf.convex-covering)

;;; `hull'
;;;
;;; Internal data structure for representing convex hulls during the
;;; decomposition computation. In contrast to the `convex-hull' type
;;; which is part of the API, this structure type
;;;
;;; 1. Stores vertex data using a concrete (double-float)
;;;    representation where `convex-hull' mirrors the component type
;;;    supplied by the caller
;;;
;;; 2. Stores additional information such as normals and bounding boxes

(defstruct (hull
            (:constructor %make-hull (vertices facets flat-p global-faces))
            (:copier NIL)
            (:predicate NIL))
  (vertices      (error "required") :type (manifolds:vertex-array manifolds:f64) :read-only T)
  (facets        (error "required") :type manifolds:face-array :read-only T)
  (facet-normals '())
  (bounding-box  NIL :type (or null cons)) ; (location . size/2) ; TODO(jmoringe): maybe store 3d-space:region directly
  ;; Maybe essential
  (flat-p        (error "required") :type boolean :read-only T)
  (global-faces  #() :type manifolds:face-array :read-only T)
  ;; Debugging
  (annotations   '() :type list)
  (problem       NIL))

(defun make-hull (vertices faces flat-p vertex-position-index)
  (let ((global-faces (make-array 0 :element-type 'manifolds:u32
                                    :adjustable T
                                    :fill-pointer 0)))
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
    (let ((global-faces (make-array (length global-faces)
                                    :element-type 'manifolds:u32
                                    :initial-contents global-faces)))
      (%make-hull vertices faces flat-p global-faces))))

(defun compute-facet-normals (hull) ; TODO(jmoringe): these are just face normals for now
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
    ;; TODO(jmoringe): 3d-spaces should accept a threshold for overlapping queries
    (center-and-size-from-min-and-max (v- min (dvec .001 .001 .001)) (v+ max (dvec .001 .001 .001)))))

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

;;; `patch'
;;;
;;; Internal data structure for representing a part of the original
;;; mesh which contributes to the convex decomposition. The patch
;;; stores the indices of mesh faces to represent its part of the
;;; original mesh. Furthermore, a valid patch stores its convex hull
;;; as a `hull' instance.
;;;
;;; The remaining slots indicate with which other patches a patch
;;; could possibly be merged as well as measures related to the cost
;;; of the merge that produced the patch.

(defstruct (patch
            (:constructor %make-patch (faces &optional surface-area hull))
            (:copier NIL)
            (:predicate NIL))
  ;; Geometry
  (faces #() :type manifolds:face-array :read-only T)
  (hull NIL :type (or null hull))
  ;; The following three slots relate to patch merging
  (links (make-array 0 :adjustable T :fill-pointer T) :type vector)
  ;; Merged surface area is always available when patches are merged,
  ;; so stored right away and never changed.
  (surface-area (error "required") :type double-float :read-only T)
  ;; Compactness is computed only if the patch is valid.
  (compactness 0.0d0 :type double-float)
  ;; Debug information
  (link))

(defun make-initial-patch (all-vertices a b c)
  (let* ((faces        (make-array 3 :element-type     'manifolds:u32
                                     :initial-contents (list a b c)))
         (surface-area (manifolds:face-area all-vertices faces 0)))
    (%make-patch faces surface-area)))

(defun compute-patch-convex-hull (all-vertices faces vertex-position-index)
  (let ((vertices (make-array 3 :element-type 'manifolds:f64
                                :adjustable T
                                :fill-pointer 0)))
    ;; Quickhull doesn't like duplicate vertices so we take care of
    ;; those here.
    (loop with seen = (make-hash-table :test #'equal)
          for j across faces
          for a of-type manifolds:u32 = (+ (* 3 j) 0) ; global vertex indices
          for b of-type manifolds:u32 = (+ (* 3 j) 1)
          for c of-type manifolds:u32 = (+ (* 3 j) 2)
          for key = (cons (logior (ash a 32) b) c)
          do (unless (gethash key seen)
               (setf (gethash key seen) T)
               (let ((x (aref all-vertices a))
                     (y (aref all-vertices b))
                     (z (aref all-vertices c)))
                 (vector-push-extend x vertices)
                 (vector-push-extend y vertices)
                 (vector-push-extend z vertices))))
    ;; Make vertex array simple.
    (let ((vertices (make-array (length vertices) :element-type 'manifolds:f64
                                                  :initial-contents vertices)))
      ;; First try to find extreme vertices than span a volume. If the
      ;; point cloud is flat, choose a different representation based
      ;; on 2D convex hull computed by the gift wrapping algorithm
      ;; (for now).
      (multiple-value-bind (base-a base-b base-c base-d)
          (let ((extrema (org.shirakumo.fraf.quickhull::compute-extrema vertices)))
            (org.shirakumo.fraf.quickhull::compute-base vertices extrema .00001d0))
        (if (null base-d)
            ;; TODO(jmoringe): do something better for the normal
            (let ((normal (vunit (manifolds:face-normal all-vertices faces 0))))
              (multiple-value-bind (vertices faces)
                  (flat-mesh vertices
                             (manifolds:v vertices base-a)
                             (manifolds:v vertices base-b)
                             (manifolds:v vertices base-c)
                             normal)
                (let ((hull (make-hull vertices faces t vertex-position-index)))
                  ;; Safety net: the normals of all hull facets should
                  ;; match the plane normal.
                  (map nil (lambda (hn)
                             (assert (< (abs (- (v. hn normal) 1)) .01)))
                       (hull-facet-normals hull))
                  hull)))
            (multiple-value-bind (vertices faces extruded-p)
                (org.shirakumo.fraf.quickhull:convex-hull
                 vertices :eps 0 :extrude-if-flat NIL)
              (assert (not extruded-p))
              (make-hull vertices faces extruded-p vertex-position-index)))))))

(defun push-link (new-link patch)
  #+assertions (assert (or (eq patch (patch-link-a new-link))
                           (eq patch (patch-link-b new-link))))
  #+assertions (assert (not (find-if (lambda (link)
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

(defun patch-debug-name (patch &key (format "~{[~d ~d ~d]~^ ~}~@[ …~]"))
  (let* ((faces (patch-faces patch))
         (count (length faces)))
    (format NIL format
            (coerce (subseq faces 0 (min 9 count)) 'list)
            (> count 9))))

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
      (let* ((hull (block nil
                     (handler-bind
                         (((or org.shirakumo.fraf.quickhull:points-colinear-error
                               org.shirakumo.fraf.quickhull:points-not-distinct-error)
                            (lambda (error)
                              (warn "Could not compute hull: ~A" error)
                              (return NIL))))
                       (compute-patch-convex-hull
                        all-vertices faces (context-vertex-position-index context)))))
             (patch (when hull (%make-patch faces surface-area hull))))
        (cond ((and patch (valid-patch-p patch context))
               (setf (patch-compactness patch) (compute-compactness all-vertices patch))
               patch)
              (T
               (values NIL patch) ; TODO(jmoringe): for debugging
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
  ;; TODO(jmoringe): When not debugging and RESULT is false, we would
  ;; not receive the patch at all, not store any result in the patch
  ;; link and not set the patch link of the patch.
  (multiple-value-bind (result debug-patch) (make-merged-patch context patch1 patch2)
    (let ((cost (when result
                  (funcall (context-cost-function context) result patch1 patch2))))
      (if cost
          (let ((link (%make-patch-link patch1 patch2 cost result)))
            (setf (patch-link result) link)
            link)
          (let ((link (%make-patch-link patch1 patch2 most-positive-double-float debug-patch)))
            (when debug-patch (setf (patch-link debug-patch) link))
            link)))))

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
    ;; Now that we actually merge this in, compute new links and
    ;; update the existing neighbour's links.
    ;;
    ;; We perform the expensive part, computation of new patches,
    ;; their hulls and the validity, in parallel. We use :ordered t to
    ;; guarantee the same result as a serial computation would
    ;; produce. This matters because the overall result depends on the
    ;; order in which patch links are updated.
    (with-tasks (:ordered t)
      (flet ((link (cur other)
               (loop with links = (patch-links cur)
                     for i from 0 below (length links)
                     for link = (aref links i)
                     do (unless (link-involves-p other link)
                          (let ((patch2 (link-other-patch link cur)))
                            (task
                             (list (make-patch-link context new-patch patch2)
                                   cur
                                   patch2)))))))
        ;; Submit tasks. Each task should compute a patch but not
        ;; touch shared data, in particular patch links.
        (link (patch-link-a link) (patch-link-b link))
        (link (patch-link-b link) (patch-link-a link))
        ;; Collect results and update patch links.
        (do-results (result)
          (destructuring-bind (new-link cur patch2) result
            (unless (find-link-involving patch2 new-patch)
              (push-link new-link new-patch))
            (let ((temp (remove-if (lambda (old-link)
                                     ;; (assert (not (link-involves-p new-patch old-link)))
                                     (when (eq old-link link)
                                       (assert (link-involves-p cur old-link)))
                                     (link-involves-p cur old-link))
                                   (patch-links patch2))))
              (setf (patch-links patch2)
                    (make-array (length temp) :initial-contents temp :adjustable T :fill-pointer T)))
            (unless (find-link-involving new-patch patch2)
              (push-link new-link patch2))))
        new-patch))))
